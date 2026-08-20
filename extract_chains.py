import os
import re
import warnings
from Bio import BiopythonWarning
from Bio.PDB import PDBParser, MMCIFParser
from Bio.SeqUtils import seq1
from Bio.SeqIO.PdbIO import PdbSeqresIterator, CifSeqresIterator
from Bio import SeqIO
from Bio.Data import PDBData
from Bio.PDB import MMCIF2Dict


def _safe_chain_id(chain_id):
    """Return (display_chain_id, filename_safe_chain_id)."""
    display = chain_id.strip() or "_"
    safe = re.sub(r"[^A-Za-z0-9_-]", "_", display)
    return display, safe


def _residue_number_string(residue):
    """
    Return the residue number as a string, combining the sequence
    identifier and insertion code if present.
    """
    # residue.id is (hetero_flag, sequence_identifier, insertion_code)
    _, seq_id, ins_code = residue.id

    base = str(seq_id)
    if ins_code and ins_code.strip():
        return base + ins_code.strip()
    return base


import os
from Bio import SeqIO
from Bio.Data import PDBData
from Bio.PDB import MMCIF2Dict


def get_seqres(path, chain_id):
    """
    Return the SEQRES sequence for the specified chain from a PDB or mmCIF file.

    Parameters
    ----------
    path : str
        Path to structure file (.pdb, .ent, .cif, .mmcif)
    chain_id : str
        Chain identifier as found in the structure file (e.g. 'A', 'B',
        auth_asym_id or label_asym_id for mmCIF).

    Returns
    -------
    str or None
        One-letter amino acid sequence of the chain, or None if not found.
    """
    ext = os.path.splitext(path)[1].lower()

    if ext in (".pdb", ".ent"):
        return _get_seqres_pdb(path, chain_id)
    elif ext in (".cif", ".mmcif", ".mcif"):
        return _get_seqres_mmcif(path, chain_id)
    else:
        raise ValueError(f"Unsupported file extension: {ext}")


def _get_seqres_pdb(path, chain_id):
    """Extract SEQRES from a PDB file using Biopython's pdb-seqres parser."""
    try:
        records = SeqIO.parse(path, "pdb-seqres")
    except Exception:
        return None

    for rec in records:
        # Biopython usually stores chain ID in annotations, or as "pdbid:chain" in id
        chain = rec.annotations.get("chain")
        if chain is None:
            # Try to extract from the record id (e.g., '1abc:A' or '1abc_A')
            if ":" in rec.id:
                chain = rec.id.split(":")[1]
            elif "_" in rec.id:
                chain = rec.id.split("_")[-1]

        if chain == chain_id:
            return str(rec.seq)

    return None


def _get_seqres_mmcif(path: str, chain_id: str):
    """Extract SEQRES from an mmCIF file using the entity_poly_seq table."""

    for record in SeqIO.parse(path, "cif-seqres"):
        if record.annotations.get("chain") == chain_id:
            return str(record.seq)


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

    res_num is the residue number from the ATOM section of the structure file
    (sequence identifier + insertion code if present).
    """
    warnings.simplefilter("ignore", BiopythonWarning)

    os.makedirs(output_dir, exist_ok=True)

    pdb_parser = PDBParser(QUIET=True)
    cif_parser = MMCIFParser(QUIET=True, auth_chains=False)

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
        best_res_nums = []  # store residue number strings
        best_res_names = []

        for chain in model:
            ca_atoms = []
            res_nums = []
            res_names = []

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
                res_nums.append(_residue_number_string(residue))
                res_names.append(seq1(residue.get_resname()))

            if len(ca_atoms) > len(best_ca_atoms):
                best_chain = chain
                best_ca_atoms = ca_atoms
                best_res_nums = res_nums
                best_res_names = res_names

        if best_chain is None or not best_ca_atoms:
            print(f"Skipping {filename}: no CA atoms found")
            continue

        chain_id_display, chain_id_safe = _safe_chain_id(best_chain.id)

        #out_filename = f"{stem}_{ext_lower.lstrip('.')}_{chain_id_safe}.txt"
        out_filename = f"{stem}.txt"
        out_path = os.path.join(output_dir, out_filename)

        with open(out_path, "w") as out_f:
            for ca, res_num, res_name in zip(best_ca_atoms, best_res_nums, best_res_names):
                x, y, z = ca.coord
                out_f.write(f"{x:.3f} {y:.3f} {z:.3f} {chain_id_display} {res_num} {res_name}\n")

        print(f"Wrote {out_filename} ({len(best_ca_atoms)} CA atoms)")

        #try:
        seq = get_seqres(path, chain_id_display)

            #seq_filename = f"{stem}_{ext_lower.lstrip('.')}_{chain_id_safe}_seq.txt"
        seq_filename = f"{stem}_seq.txt"
        seq_path = os.path.join(output_dir, seq_filename)
        seq_file = open(seq_path, 'w')
        seq_file.write(seq)
        seq_file.close()
        #except:
        #    print("Can not extract sequence of structure " + str(path) + ".")


extract_longest_chain_to_files("CIFs", "extracted_chains")
