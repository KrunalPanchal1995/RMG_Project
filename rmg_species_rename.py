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
    # keep it Cantera-friendly (no spaces)
    return name.replace(" ", "")


def rmg_from_adjacency(adj: str) -> Tuple[str, str, int, int]:
    """
    Return (smiles, formula, multiplicity, net_charge) from an RMG adjacency list.
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
    for a in mol.atoms:
        net_charge += int(getattr(a, "charge", 0))

    return smiles, str(formula), int(multiplicity), int(net_charge)


def ffcm_name_from_formula(formula: str, multiplicity: int, net_charge: int) -> Optional[str]:
    """
    Minimal mapping to match your FFCM-style list (and a few common extras).
    Extend this mapping if you want.
    """
    f = formula.replace(" ", "")

    # inerts
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
    if f == "H2O":
        return "H2O"
    if f == "HO2" and multiplicity == 2:
        return "HO2"
    if f == "H2O2":
        return "H2O2"

    # syngas
    if f == "CO" and net_charge == 0:
        return "CO"
    if f == "CO2":
        return "CO2"

    # (optional) small hydrocarbon radicals if they appear
    if f in ("CH", "CH2", "CH3", "CH4", "HCO", "CH2O"):
        return f

    return None


def uniqueify(name: str, used: set[str], suffix: str) -> str:
    """
    Ensure unique species names (required by Cantera).
    """
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


def build_mapping_from_species_dict(species_dict_path: Path, others_mode: str) -> Tuple[Dict[str, str], List[Dict[str, str]]]:
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

        # Decide new name
        ffcm = ffcm_name_from_formula(formula, mult, q)
        if ffcm is not None:
            new_name = uniqueify(ffcm, used, suffix=str(idx) if idx is not None else "x")
        else:
            if others_mode == "smiles":
                new_name = uniqueify(smiles_wo, used, suffix=str(idx) if idx is not None else "x")
            else:
                new_name = uniqueify(old_wo, used, suffix=str(idx) if idx is not None else "x")

        # IMPORTANT:
        # Only map the *exact species labels* (with number) that appear in Cantera YAML/CK,
        # to avoid touching element symbols (O, H, C) in thermo composition keys.
        mapping[old_with] = new_name

        # (Optional) map old_without_number ONLY if it cannot collide with element keys
        # Keep this conservative.
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


def replace_in_equation(equation: str, mapping: Dict[str, str]) -> str:
    # Replace only in reaction equations (safe)
    items = sorted(mapping.items(), key=lambda kv: len(kv[0]), reverse=True)
    out = equation
    for old, new in items:
        pat = r"(?<![A-Za-z0-9_])" + re.escape(old) + r"(?![A-Za-z0-9_])"
        out = re.sub(pat, new, out)
    return out


def rename_yaml_cantera(in_yaml: Path, out_yaml: Path, mapping: Dict[str, str]) -> None:
    """
    YAML-aware rename:
    - species[*]['name']
    - phases[*]['species'] list
    - reactions[*]['equation']
    - reactions[*]['efficiencies'] keys (if present)
    Does NOT touch thermo composition keys or any numeric data.
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

    # Reactions
    if "reactions" in data:
        for rxn in data["reactions"]:
            if "equation" in rxn:
                rxn["equation"] = replace_in_equation(str(rxn["equation"]), mapping)

            # third-body efficiencies keys are species names
            if "efficiencies" in rxn and isinstance(rxn["efficiencies"], dict):
                new_eff = {}
                for k, v in rxn["efficiencies"].items():
                    new_eff[mapping.get(k, k)] = v
                rxn["efficiencies"] = new_eff

    out_yaml.parent.mkdir(parents=True, exist_ok=True)
    with out_yaml.open("w", encoding="utf-8") as f:
        yaml.dump(data, f)


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species-dict", required=True, help="chemkin/species_dictionary.txt")
    ap.add_argument("--yaml", required=True, help="cantera/chem.yaml or chem_annotated.yaml")
    ap.add_argument("--out-yaml", required=True, help="output YAML path")
    ap.add_argument("--out-csv", required=True, help="species.csv path")
    ap.add_argument("--others", choices=["smiles", "keep"], default="smiles",
                    help="how to name non-FFCM species (default: smiles)")
    args = ap.parse_args()

    mapping, rows = build_mapping_from_species_dict(Path(args.species_dict), others_mode=args.others)
    rename_yaml_cantera(Path(args.yaml), Path(args.out_yaml), mapping)
    write_species_csv(rows, Path(args.out_csv))

    print(f"[OK] wrote YAML: {args.out_yaml}")
    print(f"[OK] wrote CSV:  {args.out_csv}")
    print(f"[INFO] mapped {len(rows)} species")


if __name__ == "__main__":
    main()
