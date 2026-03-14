#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from ruamel.yaml import YAML
from rmgpy.species import Species as RMGSpecies

LABEL_NUM_RE = re.compile(r"^(.*)\((\d+)\)$")


# ---------------------------------------------------------------------------
# Full FFCM1.0 target list
# ---------------------------------------------------------------------------
FFCM_TARGET_SET = {
    "AR", "HE", "N2", "H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2",
    "CO", "CO2", "C", "CH", "CH2", "CH2(S)", "CH3", "CH4", "HCO", "CH2O",
    "CH2OH", "CH3O", "CH3OH", "C2H", "C2H2", "C2H3", "C2H4", "C2H5", "C2H6",
    "HCCO", "CH2CO", "CH2CHO", "CH3CHO", "CH3CO", "H2CC", "OH*", "CH*"
}


# ---------------------------------------------------------------------------
# Manual overrides
#
# Use these for:
# - excited species like OH*, CH*
# - ambiguous isomers like CH2OH vs CH3O, CH3CO vs CH2CHO, etc.
#
# Priority:
# 1) exact RMG label with number
# 2) exact unnumbered RMG label
# 3) exact SMILES
# ---------------------------------------------------------------------------
MANUAL_OVERRIDES_BY_LABEL: Dict[str, str] = {
    # Examples:
    # "OH*(123)": "OH*",
    # "CH*(456)": "CH*",
}

MANUAL_OVERRIDES_BY_UNNUMBERED_LABEL: Dict[str, str] = {
    # Examples:
    # "OH*": "OH*",
    # "CH*": "CH*",
}

MANUAL_OVERRIDES_BY_SMILES: Dict[str, str] = {
    # Safe common examples if RMG produces exactly these smiles:
    "[CH2]O": "CH2OH",
    "C[O]": "CH3O",
    # Add more if you verify them from your species.csv output
}


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------
def split_blocks(text: str) -> List[str]:
    blocks, cur = [], []
    for line in text.splitlines():
        if line.strip() == "":
            if cur:
                blocks.append("\n".join(cur).strip("\n"))
                cur = []
        else:
            cur.append(line.rstrip("\n"))
    if cur:
        blocks.append("\n".join(cur).strip("\n"))
    return blocks


def parse_label(label: str) -> Tuple[str, str, Optional[int]]:
    label = label.strip()
    m = LABEL_NUM_RE.match(label)
    if m:
        base = m.group(1).strip()
        idx = int(m.group(2))
        return label, base, idx
    return label, label, None


def sanitize_species_name(name: str) -> str:
    # keep it cantera-friendly
    return name.replace(" ", "")


def uniqueify(name: str, used: set[str], suffix: str) -> str:
    base = sanitize_species_name(name)
    if base not in used:
        used.add(base)
        return base

    cand = f"{base}_{suffix}"
    if cand not in used:
        used.add(cand)
        return cand

    k = 2
    while True:
        cand2 = f"{base}_{suffix}_{k}"
        if cand2 not in used:
            used.add(cand2)
            return cand2
        k += 1


# ---------------------------------------------------------------------------
# RMG parsing
# ---------------------------------------------------------------------------
def rmg_from_adjacency(adj: str) -> Tuple[str, str, int, int]:
    """
    Return:
        smiles, formula, multiplicity, net_charge
    from an RMG adjacency list.
    """
    sp = RMGSpecies()
    sp.from_adjacency_list(adj)

    mol = sp.molecule[0]
    smiles = mol.to_smiles()

    try:
        formula = sp.get_formula()
    except Exception:
        formula = mol.get_formula()

    multiplicity = getattr(mol, "multiplicity", 1) or 1

    net_charge = 0
    for atom in mol.atoms:
        net_charge += int(getattr(atom, "charge", 0))

    return smiles, str(formula), int(multiplicity), int(net_charge)


# ---------------------------------------------------------------------------
# FFCM matching
# ---------------------------------------------------------------------------
def ffcm_name_from_signature(
    formula: str,
    multiplicity: int,
    net_charge: int,
    smiles: str,
    old_with: str,
    old_wo: str,
) -> Optional[str]:
    """
    Safer FFCM matcher.

    Priority:
    1) Manual override by exact RMG label
    2) Manual override by unnumbered label
    3) Manual override by exact SMILES
    4) Safe automatic mapping for unambiguous species
    """
    # -----------------
    # Manual overrides
    # -----------------
    if old_with in MANUAL_OVERRIDES_BY_LABEL:
        return MANUAL_OVERRIDES_BY_LABEL[old_with]

    if old_wo in MANUAL_OVERRIDES_BY_UNNUMBERED_LABEL:
        return MANUAL_OVERRIDES_BY_UNNUMBERED_LABEL[old_wo]

    if smiles in MANUAL_OVERRIDES_BY_SMILES:
        return MANUAL_OVERRIDES_BY_SMILES[smiles]

    # -----------------
    # Automatic mapping
    # -----------------
    f = formula.replace(" ", "")

    # Inerts
    if f == "Ar":
        return "AR"
    if f == "He":
        return "HE"
    if f == "N2":
        return "N2"

    # H/O core
    if f == "H2":
        return "H2"
    if f == "H" and multiplicity == 2:
        return "H"
    if f == "O" and multiplicity == 3:
        return "O"
    if f == "O2" and multiplicity == 3:
        return "O2"
    if f in ("HO", "OH") and multiplicity == 2:
        return "OH"
    if f == "H2O" and multiplicity == 1:
        return "H2O"
    if f == "HO2" and multiplicity == 2:
        return "HO2"
    if f == "H2O2" and multiplicity == 1:
        return "H2O2"

    # Syngas
    if f == "CO" and net_charge == 0 and multiplicity == 1:
        return "CO"
    if f == "CO2" and multiplicity == 1:
        return "CO2"

    # Carbon atom / radicals
    if f == "C" and multiplicity == 3:
        return "C"
    if f == "CH" and multiplicity == 2:
        return "CH"

    # CH2 / CH2(S)
    if f == "CH2":
        if multiplicity == 1:
            return "CH2(S)"
        if multiplicity == 3:
            return "CH2"

    if f == "CH3" and multiplicity == 2:
        return "CH3"
    if f == "CH4" and multiplicity == 1:
        return "CH4"

    if f == "CHO" and multiplicity == 2:
        return "HCO"
    if f == "CH2O" and multiplicity == 1:
        return "CH2O"

    # CH3O family (structure dependent)
    if f == "CH3O" and multiplicity == 2:
        # handled through MANUAL_OVERRIDES_BY_SMILES when possible
        return None

    if f == "CH4O" and multiplicity == 1:
        return "CH3OH"

    # C2 family
    if f == "C2H" and multiplicity == 2:
        return "C2H"
    if f == "C2H2" and multiplicity == 1:
        return "C2H2"
    if f == "C2H3" and multiplicity == 2:
        return "C2H3"
    if f == "C2H4" and multiplicity == 1:
        return "C2H4"
    if f == "C2H5" and multiplicity == 2:
        return "C2H5"
    if f == "C2H6" and multiplicity == 1:
        return "C2H6"

    # Ambiguous oxygenated C2 species:
    # HCCO, CH2CO, CH2CHO, CH3CHO, CH3CO, H2CC
    # These should be assigned via manual SMILES overrides after checking species.csv

    # Excited species OH* and CH* cannot be inferred automatically
    return None


# ---------------------------------------------------------------------------
# Mapping builder
# ---------------------------------------------------------------------------
def build_mapping_from_species_dict(
    species_dict_path: Path,
    others_mode: str
) -> Tuple[Dict[str, str], List[Dict[str, str]]]:
    """
    Build mapping old_name -> new_name, and rows for species.csv.

    others_mode:
      - 'smiles' : rename non-FFCM species to SMILES
      - 'keep'   : keep original label for non-FFCM species
    """
    text = species_dict_path.read_text(encoding="utf-8", errors="replace")
    blocks = split_blocks(text)

    used: set[str] = set()
    mapping: Dict[str, str] = {}
    rows: List[Dict[str, str]] = []

    for blk in blocks:
        lines = blk.splitlines()
        if not lines:
            continue

        label = lines[0].strip()
        old_with, old_wo, idx = parse_label(label)
        adjacency = "\n".join(lines[1:]).strip() + "\n"

        smiles_wo, formula, mult, q = rmg_from_adjacency(adjacency)

        ffcm = ffcm_name_from_signature(
            formula=formula,
            multiplicity=mult,
            net_charge=q,
            smiles=smiles_wo,
            old_with=old_with,
            old_wo=old_wo,
        )

        if ffcm is not None:
            new_name = uniqueify(ffcm, used, suffix=str(idx) if idx is not None else "x")
        else:
            if others_mode == "smiles":
                new_name = uniqueify(smiles_wo, used, suffix=str(idx) if idx is not None else "x")
            else:
                new_name = uniqueify(old_wo, used, suffix=str(idx) if idx is not None else "x")

        # IMPORTANT:
        # map exact species labels
        mapping[old_with] = new_name

        # only map unnumbered labels when safe
        if old_wo not in {"H", "O", "C", "N", "Ar", "He", "Ne"}:
            mapping.setdefault(old_wo, new_name)

        smiles_with = f"{smiles_wo}({idx})" if idx is not None else smiles_wo

        rows.append({
            "new_name": new_name,
            "rmg_label_with_number": old_with,
            "rmg_label_without_number": old_wo,
            "species_number": "" if idx is None else str(idx),
            "smiles_with_number": smiles_with,
            "smiles_without_number": smiles_wo,
            "formula": formula,
            "multiplicity": str(mult),
            "net_charge": str(q),
        })

    return mapping, rows


# ---------------------------------------------------------------------------
# Equation replacement
# ---------------------------------------------------------------------------
def replace_in_equation(equation: str, mapping: Dict[str, str]) -> str:
    items = sorted(mapping.items(), key=lambda kv: len(kv[0]), reverse=True)
    out = equation
    for old, new in items:
        pat = r"(?<![A-Za-z0-9_])" + re.escape(old) + r"(?![A-Za-z0-9_])"
        out = re.sub(pat, new, out)
    return out


# ---------------------------------------------------------------------------
# YAML-safe renaming
# ---------------------------------------------------------------------------
def rename_yaml_cantera(in_yaml: Path, out_yaml: Path, mapping: Dict[str, str]) -> None:
    """
    YAML-aware rename:
    - species[*]['name']
    - phases[*]['species'] list
    - reactions[*]['equation']
    - reactions[*]['efficiencies'] keys
    Does NOT touch thermo composition or numeric data.
    """
    yaml = YAML()
    yaml.preserve_quotes = True
    data = yaml.load(in_yaml.read_text(encoding="utf-8", errors="replace"))

    # Species blocks
    if "species" in data:
        for sp in data["species"]:
            old = sp.get("name")
            if old in mapping:
                sp["name"] = mapping[old]

    # Phase species lists
    if "phases" in data:
        for ph in data["phases"]:
            if "species" in ph and isinstance(ph["species"], list):
                ph["species"] = [mapping.get(s, s) for s in ph["species"]]

    # Reaction blocks
    if "reactions" in data:
        for rxn in data["reactions"]:
            if "equation" in rxn:
                rxn["equation"] = replace_in_equation(str(rxn["equation"]), mapping)

            if "efficiencies" in rxn and isinstance(rxn["efficiencies"], dict):
                new_eff = {}
                for k, v in rxn["efficiencies"].items():
                    new_eff[mapping.get(k, k)] = v
                rxn["efficiencies"] = new_eff

    out_yaml.parent.mkdir(parents=True, exist_ok=True)
    with out_yaml.open("w", encoding="utf-8") as f:
        yaml.dump(data, f)


# ---------------------------------------------------------------------------
# CSV writer
# ---------------------------------------------------------------------------
def write_species_csv(rows: List[Dict[str, str]], out_csv: Path) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    cols = [
        "new_name",
        "rmg_label_with_number",
        "rmg_label_without_number",
        "species_number",
        "smiles_with_number",
        "smiles_without_number",
        "formula",
        "multiplicity",
        "net_charge",
    ]
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species-dict", required=True, help="chemkin/species_dictionary.txt")
    ap.add_argument("--yaml", required=True, help="cantera/chem.yaml or chem_annotated.yaml")
    ap.add_argument("--out-yaml", required=True, help="output YAML path")
    ap.add_argument("--out-csv", required=True, help="species.csv path")
    ap.add_argument(
        "--others",
        choices=["smiles", "keep"],
        default="smiles",
        help="how to name non-FFCM species (default: smiles)"
    )
    args = ap.parse_args()

    mapping, rows = build_mapping_from_species_dict(
        Path(args.species_dict),
        others_mode=args.others
    )
    rename_yaml_cantera(Path(args.yaml), Path(args.out_yaml), mapping)
    write_species_csv(rows, Path(args.out_csv))

    print(f"[OK] wrote YAML: {args.out_yaml}")
    print(f"[OK] wrote CSV:  {args.out_csv}")
    print(f"[INFO] mapped {len(rows)} species")


if __name__ == "__main__":
    main()
