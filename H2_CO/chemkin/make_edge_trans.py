#!/usr/bin/env python3
import argparse
import os
from pathlib import Path

from rmgpy import settings
from rmgpy.chemkin import load_species_dictionary, save_transport_file
from rmgpy.data.transport import TransportDatabase
import rmgpy.data.rmg as rmg_data


def ensure_global_transport_db(tdb: TransportDatabase):
    """
    RMG's Species.get_transport_data() looks up the transport DB via rmgpy.data.rmg.get_db('transport').
    That requires rmgpy.data.rmg.database.transport to exist.
    """
    if getattr(rmg_data, "database", None) is None:
        class DummyDB:
            pass
        rmg_data.database = DummyDB()

    rmg_data.database.transport = tdb


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species-dict", required=True, help="Path to species_edge_dictionary.txt")
    ap.add_argument("--out", required=True, help="Output transport file path (e.g., tran_edge.dat)")
    ap.add_argument(
        "--transport-library",
        default=None,
        help="Transport library label (case-sensitive), e.g. PrimaryTransportLibrary",
    )
    args = ap.parse_args()

    species_dict_path = Path(args.species_dict)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # 1) Load edge species (core+edge species list)
    sp_dict = load_species_dictionary(str(species_dict_path), generate_resonance_structures=True)
    species_list = list(sp_dict.values())
    print(f"[INFO] Loaded {len(species_list)} species from {species_dict_path}")

    # 2) Load transport database
    db_dir = settings["database.directory"]
    transport_path = os.path.join(db_dir, "transport")

    tdb = TransportDatabase()
    if args.transport_library:
        # library name is case-sensitive
        tdb.load(transport_path, libraries=[args.transport_library])
        print(f"[INFO] Loaded transport library: {args.transport_library}")
    else:
        # load all transport libraries if not specified
        tdb.load(transport_path)
        print("[INFO] Loaded all available transport libraries")

    # 3) Register the transport DB globally so Species.get_transport_data() can find it
    ensure_global_transport_db(tdb)

    # 4) Write transport file for the full edge species list
    save_transport_file(str(out_path), species_list)
    print(f"[OK] Wrote tran file: {out_path} (species={len(species_list)})")


if __name__ == "__main__":
    main()
