import os
import re
import warnings
from Bio import BiopythonWarning
from Bio.PDB import MMCIFParser
from Bio.SeqUtils import seq1
from Bio.SeqIO.PdbIO import CifSeqresIterator
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


def get_seqres(path, chain_id):
    """
    Return the SEQRES sequence for the specified chain from a mmCIF file.

    Parameters
    ----------
    path : str
        Path to structure file.
    chain_id : str
        Chain identifier as found in the structure file.

    Returns
    -------
    str or None
        One-letter amino acid sequence of the chain, or None if not found.
    """

    ext = os.path.splitext(path)[1].lower()

    if ext not in (".cif", ".mmcif", ".mcif"):
        raise ValueError(f"Only CIF files are supported: {path}")

    return _get_seqres_mmcif(path, chain_id)


def _get_seqres_mmcif(path: str, chain_id: str):
    """
    Extract the SEQRES-like sequence from an mmCIF file.

    If multiple residues have the same sequence position,
    only the first encountered residue is kept.
    """

    cif = MMCIF2Dict.MMCIF2Dict(path)

    try:
        entity_ids = cif["_entity_poly_seq.entity_id"]
        num = cif["_entity_poly_seq.num"]
        mon_id = cif["_entity_poly_seq.mon_id"]
    except KeyError:
        return None

    # Mapping from entity_id to chain IDs
    try:
        poly_entity_ids = cif["_struct_asym.entity_id"]
        poly_chain_ids = cif["_struct_asym.id"]
    except KeyError:
        return None

    entity_for_chain = None

    for entity_id, asym_id in zip(poly_entity_ids, poly_chain_ids):
        if asym_id == chain_id:
            entity_for_chain = entity_id
            break

    if entity_for_chain is None:
        return None

    sequence = []
    seen_positions = set()

    for entity_id, position, residue_name in zip(
        entity_ids,
        num,
        mon_id
    ):
        if entity_id != entity_for_chain:
            continue

        if position in seen_positions:
            continue

        seen_positions.add(position)

        sequence.append(seq1(residue_name))

    return "".join(sequence)


def get_first_residue_variants(path, chain_id):
    """
    Return the first residue variant encountered in the mmCIF file
    for each residue position of the specified chain.

    If several residue variants occur at the same position,
    only the first one encountered in the file is kept.

    Parameters
    ----------
    path : str
        Path to the mmCIF file.
    chain_id : str
        Chain identifier.

    Returns
    -------
    dict
        Mapping from residue number to the first residue name
        encountered at that position.
    """

    cif = MMCIF2Dict.MMCIF2Dict(path)

    try:
        asym_ids = cif["_atom_site.label_asym_id"]
        seq_ids = cif["_atom_site.label_seq_id"]
        comp_ids = cif["_atom_site.label_comp_id"]
    except KeyError:
        return {}

    ins_codes = cif.get(
        "_atom_site.pdbx_PDB_ins_code",
        ["."] * len(asym_ids)
    )

    first_variants = {}

    for asym_id, seq_id, ins_code, comp_id in zip(
        asym_ids,
        seq_ids,
        ins_codes,
        comp_ids
    ):
        if asym_id != chain_id:
            continue

        if seq_id in (".", "?"):
            continue

        if ins_code in (".", "?"):
            ins_code = ""

        residue_number = str(seq_id) + str(ins_code).strip()

        # Keep only the FIRST residue encountered
        # at this position in the original mmCIF file.
        if residue_number not in first_variants:
            first_variants[residue_number] = comp_id

    return first_variants


def extract_longest_chain_to_files(input_dir, output_dir):
    """
    Process all mmCIF structures in input_dir.

    For each structure:
      - parse the first model
      - select the chain with the most CA atoms
      - if a residue position has multiple residue variants,
        select the variant that appears first in the original
        mmCIF file
      - write only the selected CA atoms to output_dir
      - extract and write the SEQRES-like sequence

    Output filename:
      <pdbid>.txt

    Output line format:
      x y z chain_id res_num res_name
    """

    warnings.simplefilter("ignore", BiopythonWarning)

    os.makedirs(output_dir, exist_ok=True)

    cif_parser = MMCIFParser(QUIET=True, auth_chains=False)
    cif_exts = {".cif", ".mmcif", ".mcif"}

    for filename in sorted(os.listdir(input_dir)):

        path = os.path.join(input_dir, filename)

        if not os.path.isfile(path):
            continue

        stem, ext = os.path.splitext(filename)
        ext_lower = ext.lower()

        if ext_lower not in cif_exts:
            continue

        try:
            structure = cif_parser.get_structure(stem, path)
        except Exception as exc:
            print(f"Skipping {filename}: {exc}")
            continue

        if len(structure) == 0:
            print(f"Skipping {filename}: no models")
            continue

        model = structure[0]

        best_chain = None
        best_ca_count = 0

        for chain in model:

            ca_count = 0

            for residue in chain:

                if "CA" in residue:
                    ca_count += 1

            if ca_count > best_ca_count:
                best_chain = chain
                best_ca_count = ca_count

        if best_chain is None or best_ca_count == 0:
            print(f"Skipping {filename}: no CA atoms found")
            continue

        chain_id_display, chain_id_safe = _safe_chain_id(best_chain.id)

        first_variants = get_first_residue_variants(
            path,
            chain_id_display
        )

        best_ca_atoms = []
        best_res_nums = []
        best_res_names = []

        for residue in best_chain:

            if "CA" not in residue:
                continue

            res_num = _residue_number_string(residue)

            # Which residue appeared FIRST in the original mmCIF?
            first_res_name = first_variants.get(
                res_num,
                residue.get_resname()
            )

            selected_residue = None

            if not residue.is_disordered():

                if residue.get_resname() == first_res_name:
                    selected_residue = residue

            # Disordered residue:
            # several amino acids at the same position
            else:

                try:
                    for variant in residue.disordered_get_list():

                        if variant.get_resname() == first_res_name:
                            selected_residue = variant
                            break

                except AttributeError:
                    pass

            if selected_residue is None:

                if not residue.is_disordered():
                    selected_residue = residue

                else:
                    try:
                        selected_residue = residue.disordered_get_list()[0]
                    except (AttributeError, IndexError):
                        continue

            if "CA" not in selected_residue:
                continue

            ca = selected_residue["CA"]

            # CA itself can have alternative locations.
            # Here we take the first CA object belonging
            # to the already selected residue.
            if ca.is_disordered():

                try:
                    ca = ca.disordered_get_list()[0]
                except (AttributeError, IndexError):
                    continue

            best_ca_atoms.append(ca)
            best_res_nums.append(res_num)
            best_res_names.append(
                seq1(selected_residue.get_resname())
            )


        if not best_ca_atoms:
            print(f"Skipping {filename}: no selected CA atoms")
            continue

        out_filename = f"{stem}.txt"
        out_path = os.path.join(output_dir, out_filename)

        with open(out_path, "w") as out_f:

            for ca, res_num, res_name in zip(
                best_ca_atoms,
                best_res_nums,
                best_res_names
            ):
                x, y, z = ca.coord

                out_f.write(
                    f"{x:.3f} {y:.3f} {z:.3f} "
                    f"{chain_id_display} {res_num} {res_name}\n"
                )

        print(
            f"Wrote {out_filename} "
            f"({len(best_ca_atoms)} CA atoms)"
        )

        seq = get_seqres(path, chain_id_display)

        if seq is None:
            print(
                f"SEQRES not found for chain "
                f"{chain_id_display} in {filename}"
            )
            continue

        seq_filename = f"{stem}_seq.txt"
        seq_path = os.path.join(output_dir, seq_filename)

        with open(seq_path, "w") as seq_file:
            seq_file.write(seq)


extract_longest_chain_to_files("CIFs", "extracted_chains")
