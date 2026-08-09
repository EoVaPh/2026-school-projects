from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.vectors import calc_dihedral
from pathlib import Path
from Bio.PDB.vectors import calc_dihedral
import traceback


amino_acids = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
}


def get_amino_acid(code):
    return amino_acids.get(code.upper(), "X")


def get_phi_psi(filename: Path, output_file: Path):
    """
    Calculate phi and psi dihedral angles for residues in the first
    chain of a PDB or mmCIF structure and append the results to a
    common output file.

    Alternative residue locations are handled by selecting the first
    available residue variant. Alternative atom locations are also
    handled by selecting the first available atom variant.

    Non-canonical amino acids are represented as X.
    Residues without a CA atom are skipped.

    Parameters:
        filename (Path): Path to the input PDB or mmCIF file.
        output_file (Path): Path to the common output file.

    The angles are stored in radians.
    """

    filename = Path(filename)
    output_file = Path(output_file)

    if filename.suffix.lower() == ".pdb":
        parser = PDBParser(QUIET=True)
    elif filename.suffix.lower() in [".cif", ".mmcif"]:
        parser = MMCIFParser(QUIET=True)
    else:
        raise ValueError(f"Unsupported file format: {filename}")

    structure = parser.get_structure("structure", filename)
    model = next(structure.get_models())
    chain = next(model.get_chains())

    print(f"Processing: {filename.name}")
    print(f"First chain: {chain.id}")


    def get_residue_variant(residue):

        if residue.is_disordered() == 2:
            variants = list(residue.child_dict.values())
            if variants:
                return variants[0]
            return None
        
        return residue


    def get_atom(residue, atom_name):
        if residue is None:
            return None
        if atom_name not in residue:
            return None
        atom = residue[atom_name]

        if atom.is_disordered():
            variants = list(atom.child_dict.values())
            if variants:
                return variants[0]
            return None

        return atom


    residues = []

    for residue in chain.get_residues():
        selected_residue = get_residue_variant(residue)

        if selected_residue is None:
            continue
        if get_atom(selected_residue, "CA") is None:
            continue

        residues.append(selected_residue)

    with open(output_file, "a", encoding="utf-8") as file:

        file.write(f">{filename.name}\n")

        for i, residue in enumerate(residues):
            amino_acid = get_amino_acid(residue.resname)

            residue_number = residue.id[1]

            phi = None
            psi = None

            if i > 0:
                previous_residue = residues[i - 1]
                previous_C = get_atom(previous_residue, "C")
                current_N = get_atom(residue, "N")
                current_CA = get_atom(residue, "CA")
                current_C = get_atom(residue, "C")

                if all([previous_C is not None, current_N is not None, current_CA is not None, current_C is not None]):

                    phi = calc_dihedral(previous_C.get_vector(), current_N.get_vector(), current_CA.get_vector(), current_C.get_vector())

            if i < len(residues) - 1:
                next_residue = residues[i + 1]
                current_N = get_atom(residue, "N")
                current_CA = get_atom(residue, "CA")
                current_C = get_atom(residue, "C")
                next_N = get_atom(next_residue, "N")

                if all([current_N is not None, current_CA is not None, current_C is not None, next_N is not None]):

                    psi = calc_dihedral(current_N.get_vector(), current_CA.get_vector(), current_C.get_vector(), next_N.get_vector())

            file.write(f"{chain.id}\t"
                f"{residue_number}\t"
                f"{amino_acid}\t"
                f"{phi}\t"
                f"{psi}\n"
            )

        file.write("\n")


def check_sequence(angles_file: Path, structure_name: str, sequence_file: Path) -> bool:

    angles_seq = ""
    inside_structure = False

    with open(angles_file, "r", encoding="utf-8") as file:
        for line in file:
            line = line.strip()

            if not line:
                continue
            if line.startswith(">"):
                current_name = line[1:].strip()
                inside_structure = (current_name == structure_name)
                continue

            if not inside_structure:
                continue

            parts = line.split()

            if len(parts) < 3:
                continue

            angles_seq += parts[2]

    atom_seq = ""

    with open(sequence_file, "r", encoding="utf-8") as file:

        structure_lines = file.readlines()
        atom_seq = structure_lines[3]

    angles_seq = angles_seq.strip().upper()
    atom_seq = atom_seq.strip().upper()

    if angles_seq != atom_seq:

        print(f"\nSequence mismatch: {structure_name}")
        print(f"angles_seq: {angles_seq}")
        print(f"atom_seq:   {atom_seq}")
        print(f"angles length: {len(angles_seq)}")
        print(f"atom length:   {len(atom_seq)}")

        min_length = min(len(angles_seq), len(atom_seq))

        for i in range(min_length):
            if angles_seq[i] != atom_seq[i]:
                print(
                    f"First difference at position "
                    f"{i + 1}: "
                    f"{angles_seq[i]} != {atom_seq[i]}"
                )
                break

        return False

    return True


structure_folder = Path("amyloid_structures")
sequences_folder = Path("sequences")
output_file = Path("all_phi_psi.txt")

output_file.write_text("", encoding="utf-8")

for structure_file in sorted(structure_folder.iterdir()):

    if (not structure_file.is_file() or structure_file.suffix.lower() not in [".pdb", ".cif", ".mmcif"]):
        continue

    try:
        get_phi_psi(structure_file, output_file)
        sequence_file = (sequences_folder/ f"{structure_file.name}_sequences.fa")

        if sequence_file.exists():

            result = check_sequence(output_file, structure_file.name, sequence_file)

            if result:
                print(
                    f"Sequence check True for "
                    f"{structure_file}"
                )

            else:
                print(
                    f"Sequence check False for "
                    f"{structure_file}"
                )

        else:

            print(
                f"Sequence file not found: "
                f"{sequence_file}"
            )

    except Exception as e:

        print(
            f"ERROR processing "
            f"{structure_file.name}: {e}"
        )

        traceback.print_exc()


print("Completed")
