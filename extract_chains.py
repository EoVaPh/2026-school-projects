import os
import re
import warnings
from Bio import BiopythonWarning
from Bio.PDB import PDBParser, MMCIFParser


def _safe_chain_id(chain_id):
    """Return (display_chain_id, filename_safe_chain_id)."""
    display = chain_id.strip() or "_"
    safe = re.sub(r"[^A-Za-z0-9_-]", "_", display)
    return display, safe


def extract_longest_chain_to_files(input_dir, output_dir):
    """
    Process all PDB/mmCIF structures in input_dir.

    For each structure:
      - parse the first model
      - select the chain with the most CA atoms
      - write CA coordinates to output_dir

    Output filename:
      <pdbid>_<extension>_<chain>.txt

    Output line format:
      x y z chain_id res_num

    res_num is the 1-based sequential position of the residue in the chain.
    """
    warnings.simplefilter("ignore", BiopythonWarning)

    os.makedirs(output_dir, exist_ok=True)

    pdb_parser = PDBParser(QUIET=True)
    cif_parser = MMCIFParser(QUIET=True)

    pdb_exts = {".pdb", ".ent"}
    cif_exts = {".cif", ".mmcif", ".mcif"}

    for filename in sorted(os.listdir(input_dir)):
        path = os.path.join(input_dir, filename)
        if not os.path.isfile(path):
            continue

        stem, ext = os.path.splitext(filename)
        ext_lower = ext.lower()

        if ext_lower in pdb_exts:
            parser = pdb_parser
        elif ext_lower in cif_exts:
            parser = cif_parser
        else:
            continue

        try:
            structure = parser.get_structure(stem, path)
        except Exception as exc:
            print(f"Skipping {filename}: {exc}")
            continue

        if len(structure) == 0:
            print(f"Skipping {filename}: no models")
            continue

        model = structure[0]

        best_chain = None
        best_ca_atoms = []

        for chain in model:
            ca_atoms = []

            for residue in chain:
                if "CA" not in residue:
                    continue

                ca = residue["CA"]

                # If the CA atom has alternate locations, take the first one.
                try:
                    if ca.is_disordered():
                        ca = ca.disordered_get_list()[0]
                except AttributeError:
                    pass

                ca_atoms.append(ca)

            if len(ca_atoms) > len(best_ca_atoms):
                best_chain = chain
                best_ca_atoms = ca_atoms

        if best_chain is None or not best_ca_atoms:
            print(f"Skipping {filename}: no CA atoms found")
            continue

        chain_id_display, chain_id_safe = _safe_chain_id(best_chain.id)

        out_filename = f"{stem}_{ext_lower.lstrip('.')}_{chain_id_safe}.txt"
        out_path = os.path.join(output_dir, out_filename)

        with open(out_path, "w") as out_f:
            for idx, ca in enumerate(best_ca_atoms, start=1):
                x, y, z = ca.coord
                out_f.write(f"{x:.3f} {y:.3f} {z:.3f} {chain_id_display} {idx}\n")

        print(f"Wrote {out_filename} ({len(best_ca_atoms)} CA atoms)")

extract_longest_chain_to_files("amyloid_structures", "extracted_chains")
