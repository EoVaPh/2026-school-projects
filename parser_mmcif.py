from pathlib import Path
from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


amino_acids = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
}


def get_amino_acid(code):
    return amino_acids.get(code.upper(), "X")


folder = Path(
    r"C:\Users\User\Documents\bioinf\smtb\all_amyloid_structures\amyloid_structures"
)

parser = MMCIFParser(QUIET=True)


for structure_file in list(folder.glob("*.cif")) + list(folder.glob("*.mmcif")):
    name = structure_file.stem
    print(f"Processing: {structure_file.name}")

    try:
        structure = parser.get_structure(name, structure_file)
        mmcif = MMCIF2Dict(str(structure_file))
    except Exception as error:
        print(f"Error reading {structure_file.name}: {error}")
        continue

    records = []
    seq = ""
    seq_res = ""

    # Use the first model
    model = next(structure.get_models(), None)

    if model is None:
        print(f"No models found: {structure_file.name}")
        continue

    # Use the first chain
    chain = next(model.get_chains(), None)

    if chain is None:
        print(f"No chains found: {structure_file.name}")
        continue

    chain_id = chain.id

    # Get the entity ID corresponding to the first chain
    seq_res = ""

    entity_ids = mmcif.get("_struct_asym.entity_id", [])
    seq_entity_ids = mmcif.get("_entity_poly_seq.entity_id", [])
    seq_mon_ids = mmcif.get("_entity_poly_seq.mon_id", [])
    seq_nums = mmcif.get("_entity_poly_seq.num", [])

    if entity_ids and seq_entity_ids:
        entity_id = entity_ids[0]

        sequence_data = []

        for current_entity_id, num, mon_id in zip(
            seq_entity_ids, seq_nums, seq_mon_ids
        ):
            if current_entity_id == entity_id:
                sequence_data.append((int(num), mon_id))

        sequence_data.sort()

        for _, mon_id in sequence_data:
            seq_res += get_amino_acid(mon_id)
            
    # Get the complete sequence from _entity_poly_seq
    if entity_id is not None:
        seq_entity_ids = mmcif.get("_entity_poly_seq.entity_id", [])
        seq_nums = mmcif.get("_entity_poly_seq.num", [])
        seq_mon_ids = mmcif.get("_entity_poly_seq.mon_id", [])

        if isinstance(seq_entity_ids, str):
            seq_entity_ids = [seq_entity_ids]

        if isinstance(seq_nums, str):
            seq_nums = [seq_nums]

        if isinstance(seq_mon_ids, str):
            seq_mon_ids = [seq_mon_ids]

        sequence_data = []

        for current_entity_id, num, mon_id in zip(
            seq_entity_ids,
            seq_nums,
            seq_mon_ids
        ):
            if current_entity_id == entity_id:
                try:
                    sequence_data.append((int(num), mon_id))
                except ValueError:
                    continue

        sequence_data.sort(key=lambda x: x[0])

        for _, mon_id in sequence_data:
            seq_res += get_amino_acid(mon_id)

    # Get CA atoms from the first chain
    for residue in chain:
        resname = residue.get_resname().strip()
        amino_acid = get_amino_acid(resname)

        # Skip residues without a CA atom
        if "CA" not in residue:
            continue

        ca = residue["CA"]

        # Check alternate location
        altloc = ca.get_altloc()

        if altloc not in (" ", "", "A"):
            continue

        # Get CA coordinates
        x, y, z = ca.get_coord()

        # Get residue number
        residue_number = residue.id[1]

        one_record = (
            float(x),
            float(y),
            float(z),
            amino_acid,
            int(residue_number)
        )

        records.append(one_record)
        seq += amino_acid

    # Write sequences
    with open(folder.parent / f"{name}_sequences.fa", "w") as file:
        print(">SEQRES", file=file)
        print(seq_res, file=file)
        print(">ATOM sequence", file=file)
        print(seq, file=file)

    # Write CA coordinates
    with open(folder.parent / f"{name}_records.txt", "w") as file:
        for record in records:
            for item in record:
                print(item, end="\t", file=file)
            print(file=file)
