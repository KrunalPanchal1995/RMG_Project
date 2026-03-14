#!/usr/bin/env python3
"""
Rename RMG-generated Cantera YAML species to more human-readable names using the RMG species dictionary.

Inputs
------
1) species dictionary text file (RMG style), e.g. chemkin/species_dictionary.txt
2) cantera YAML file, e.g. cantera/chem.yaml or cantera/chem_annotated.yaml

Outputs
-------
1) New Cantera YAML with renamed species (all occurrences updated consistently)
2) species.csv cross-reference:
   - new_name
   - old_with_number
   - old_without_number
   - species_number
   - formula
   - multiplicity
   - net_charge

How "human-readable" names are made
-----------------------------------
- We compute the molecular formula from the adjacency list (Hill order).
- Base new name = formula (e.g., HO2, CHO2, C2H2O4, etc.)
- If multiple species share the same formula, we disambiguate using multiplicity and/or the RMG index:
    formula_m{mult} or formula_m{mult}_{id}

This ensures names are:
- readable,
- deterministic,
- unique (required by Cantera).

Usage
-----
python rename_rmg_cantera_species.py \
    --species-dict chemkin/species_dictionary.txt \
    --yaml cantera/chem.yaml \
    --out-yaml cantera/chem_readable.yaml \
    --out-csv cantera/species.csv
"""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple


@dataclass
class SpeciesEntry:
    old_with_number: str
    old_without_number: str
    species_number: int | None
    formula: str
    multiplicity: int
    net_charge: int
    new_name: str = ""


LABEL_NUM_RE = re.compile(r"^(.*)\((\d+)\)$")


def split_blocks(text: str) -> List[str]:
    # Species dictionary blocks are separated by blank lines
    blocks = []
    cur: List[str] = []
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


def parse_label(label: str) -> Tuple[str, str, int | None]:
    """Return (old_with_number, old_without_number, species_number)."""
    label = label.strip()
    m = LABEL_NUM_RE.match(label)
    if m:
        base = m.group(1).strip()
        idx = int(m.group(2))
        return label, base, idx
    return label, label, None


def parse_multiplicity_and_atoms(block_lines: List[str]) -> Tuple[int, Dict[str, int], int]:
    """
    Parse multiplicity, element counts, and net charge from the adjacency list.
    Lines look like:
        multiplicity 3
        1 O u1 p2 c0 {2,S}
    """
    multiplicity = 1
    counts: Dict[str, int] = {}
    net_charge = 0

    for ln in block_lines[1:]:
        s = ln.strip()
        if not s:
            continue
        if s.lower().startswith("multiplicity"):
            parts = s.split()
            if len(parts) >= 2:
                try:
                    multiplicity = int(parts[1])
                except ValueError:
                    pass
            continue

        # Atom line: starts with an integer index
        parts = s.split()
        if not parts:
            continue
        if not parts[0].isdigit():
            continue

        # element is the second token (e.g., O, C, H, Ar, He)
        if len(parts) < 2:
            continue
        elem = parts[1]
        counts[elem] = counts.get(elem, 0) + 1

        # find a token like c0, c+1, c-1
        for tok in parts:
            if tok.startswith("c") and len(tok) >= 2:
                # tok can be c0, c+1, c-1
                try:
                    net_charge += int(tok[1:])
                except ValueError:
                    pass
                break

    return multiplicity, counts, net_charge


def hill_formula(counts: Dict[str, int]) -> str:
    """
    Hill system: C, H, then alphabetical (Ar, He, N, O, ...).
    """
    def fmt(elem: str, n: int) -> str:
        return f"{elem}{n}" if n != 1 else elem

    c = counts.get("C", 0)
    h = counts.get("H", 0)

    out = ""
    used = set()

    if c > 0:
        out += fmt("C", c)
        used.add("C")
    if h > 0:
        out += fmt("H", h)
        used.add("H")

    # remaining elements in alphabetical order
    for elem in sorted(k for k in counts.keys() if k not in used):
        out += fmt(elem, counts[elem])

    # If somehow empty (shouldn't happen), fallback
    return out or "X"


def build_entries(species_dict_path: Path) -> List[SpeciesEntry]:
    text = species_dict_path.read_text(encoding="utf-8", errors="replace")
    blocks = split_blocks(text)

    entries: List[SpeciesEntry] = []
    for blk in blocks:
        lines = blk.splitlines()
        if not lines:
            continue

        label = lines[0].strip()
        old_with, old_wo, idx = parse_label(label)

        mult, counts, q = parse_multiplicity_and_atoms(lines)
        formula = hill_formula(counts)

        entries.append(
            SpeciesEntry(
                old_with_number=old_with,
                old_without_number=old_wo,
                species_number=idx,
                formula=formula,
                multiplicity=mult,
                net_charge=q,
            )
        )
    return entries


def assign_unique_new_names(entries: List[SpeciesEntry]) -> None:
    """
    Base name = formula.
    If duplicates exist, disambiguate by multiplicity and/or species number.
    """
    # group by formula
    by_formula: Dict[str, List[SpeciesEntry]] = {}
    for e in entries:
        by_formula.setdefault(e.formula, []).append(e)

    used_names: set[str] = set()

    for formula, group in by_formula.items():
        if len(group) == 1:
            name = formula
            # ensure uniqueness globally
            if name in used_names:
                # extremely rare: fallback
                suffix = group[0].species_number if group[0].species_number is not None else id(group[0])
                name = f"{formula}_{suffix}"
            group[0].new_name = name
            used_names.add(name)
            continue

        # For duplicates: try formula_m{mult} first
        # then add _{id} if still duplicates.
        temp_names: Dict[str, int] = {}
        for e in group:
            base = f"{formula}_m{e.multiplicity}"
            temp_names[base] = temp_names.get(base, 0) + 1

        for e in group:
            base = f"{formula}_m{e.multiplicity}"
            if temp_names[base] == 1 and base not in used_names:
                e.new_name = base
            else:
                # add species number for guaranteed uniqueness
                sid = e.species_number if e.species_number is not None else "x"
                cand = f"{base}_{sid}"
                # final fallback if somehow still collides
                k = 2
                final = cand
                while final in used_names:
                    final = f"{cand}_{k}"
                    k += 1
                e.new_name = final
            used_names.add(e.new_name)


def build_mapping(entries: List[SpeciesEntry]) -> Dict[str, str]:
    """
    Map *exact* old labels in YAML to new names.
    We map the 'old_with_number' form, because that is what RMG typically writes to YAML.
    Also map old_without_number if it differs (helps some cases).
    """
    m: Dict[str, str] = {}
    for e in entries:
        m[e.old_with_number] = e.new_name
        if e.old_without_number != e.old_with_number:
            # only add if not already mapped
            m.setdefault(e.old_without_number, e.new_name)
    return m


def safe_token_replace(text: str, mapping: Dict[str, str]) -> str:
    """
    Replace species names using regex token boundaries to avoid partial replacements,
    e.g. avoid changing CO2 when replacing CO.

    Boundary rule:
      old must not be preceded by [A-Za-z0-9_]
      old must not be followed by   [A-Za-z0-9_]
    This works well for Cantera YAML, where species are separated by spaces, commas, +, <=>, :, quotes, brackets, etc.
    """
    # Replace longer names first
    items = sorted(mapping.items(), key=lambda kv: len(kv[0]), reverse=True)

    out = text
    for old, new in items:
        pattern = r"(?<![A-Za-z0-9_])" + re.escape(old) + r"(?![A-Za-z0-9_])"
        out = re.sub(pattern, new, out)
    return out


def write_species_csv(entries: List[SpeciesEntry], out_csv: Path) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([
            "new_name",
            "old_with_number",
            "old_without_number",
            "species_number",
            "formula",
            "multiplicity",
            "net_charge",
        ])
        for e in sorted(entries, key=lambda x: (x.formula, x.species_number or -1, x.old_with_number)):
            w.writerow([
                e.new_name,
                e.old_with_number,
                e.old_without_number,
                "" if e.species_number is None else e.species_number,
                e.formula,
                e.multiplicity,
                e.net_charge,
            ])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--species-dict", required=True, help="Path to RMG species_dictionary.txt")
    ap.add_argument("--yaml", required=True, help="Path to Cantera YAML (chem.yaml or chem_annotated.yaml)")
    ap.add_argument("--out-yaml", required=True, help="Output Cantera YAML path with renamed species")
    ap.add_argument("--out-csv", required=True, help="Output species.csv cross-reference path")
    args = ap.parse_args()

    species_dict_path = Path(args.species_dict)
    yaml_path = Path(args.yaml)
    out_yaml = Path(args.out_yaml)
    out_csv = Path(args.out_csv)

    entries = build_entries(species_dict_path)
    if not entries:
        raise RuntimeError(f"No species parsed from {species_dict_path}")

    assign_unique_new_names(entries)
    mapping = build_mapping(entries)

    yaml_text = yaml_path.read_text(encoding="utf-8", errors="replace")
    new_yaml_text = safe_token_replace(yaml_text, mapping)

    out_yaml.parent.mkdir(parents=True, exist_ok=True)
    out_yaml.write_text(new_yaml_text, encoding="utf-8")

    write_species_csv(entries, out_csv)

    print(f"[OK] Wrote renamed Cantera YAML: {out_yaml}")
    print(f"[OK] Wrote cross-reference CSV: {out_csv}")
    print(f"[INFO] Species renamed: {len(entries)}")


if __name__ == "__main__":
    main()
