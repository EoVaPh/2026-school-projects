import os
from pathlib import Path
from typing import Optional
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser


def longest_chain_id(structure_path: str) -> Optional[str]:
    """
    Return the chain ID of the longest chain in a PDB or CIF structure file.

    Chain length is measured as the number of ATOM residues (polymer residues)
    in the first model. HETATM records (ligands, water, etc.) are ignored.

    Requires Biopython:  pip install biopython
    """
    path = Path(structure_path)
    suffix = path.suffix.lower()

    if suffix == ".pdb":
        parser = PDBParser(QUIET=True)
    elif suffix in {".cif", ".mmcif"}:
        parser = MMCIFParser(QUIET=True)
    else:
        raise ValueError("Unsupported file extension. Use .pdb or .cif")

    structure = parser.get_structure("structure", str(path))
    model = structure[0]

    best_chain_id = None
    best_length = -1

    for chain in model:
        # Count only standard polymer residues (ATOM records)
        n_residues = sum(1 for res in chain if res.id[0] == " ")

        if n_residues > best_length:
            best_length = n_residues
            best_chain_id = chain.id

    return best_chain_id


def pre_calculate(pdb_ids: list):
    pdb_structures = []

    structures_path = Path('../amyloid_structures')

    file_names = set([f.name for f in structures_path.iterdir() if f.is_file()])

    for pdb_id in pdb_ids:
        if pdb_id + '.cif' in file_names:
            pdb_structures.append(pdb_id + '.cif')
        else:
            pdb_structures.append(pdb_id + '.pdb')

    counter = 0

    for structure in pdb_structures:
        output = 'chain_' + structure

        chain = longest_chain_id('../amyloid_structures/' + structure)

        helix_command_template = 'python helix.py STRUCTURE --chains CHAIN --twist 0.00 --rise 4.75 --box_and_pixel 384 1.067 --n_copies 4 --output OUTPUT'

        os.system(helix_command_template.replace('STRUCTURE', '../amyloid_structures/' + structure).replace('OUTPUT', output).replace('CHAIN', chain))

        os.system('python contacts.py OUTPUT'.replace('OUTPUT', output))

        contacts = output[:-4] + '_contacts.csv'

        counter += 1

        print(counter, '/', len(pdb_structures))


def families_apd(pdb_ids: list, family_path: str):
    pdb_structures = []

    structures_path = Path('../amyloid_structures')

    file_names = set([f.name for f in structures_path.iterdir() if f.is_file()])

    for pdb_id in pdb_ids:
        if pdb_id + '.cif' in file_names:
            pdb_structures.append(pdb_id + '.cif')
        else:
            pdb_structures.append(pdb_id + '.pdb')

    for structure_1 in pdb_structures:
        for structure_2 in pdb_structures:
            if structure_1 >= structure_2: continue

            output_1 = 'chain_' + structure_1
            output_2 = 'chain_' + structure_2

            #chain_1 = longest_chain_id('../amyloid_structures/' + structure_1)
            #chain_2 = longest_chain_id('../amyloid_structures/' + structure_2)

            #helix_command_template = 'python helix.py STRUCTURE --chains CHAIN --twist 0.00 --rise 4.75 --box_and_pixel 384 1.067 --n_copies 4 --output OUTPUT'

            #os.system(helix_command_template.replace('STRUCTURE', '../amyloid_structures/' + structure_1).replace('OUTPUT', output_1).replace('CHAIN', chain_1))
            #os.system(helix_command_template.replace('STRUCTURE', '../amyloid_structures/' + structure_2).replace('OUTPUT', output_2).replace('CHAIN', chain_2))

            #os.system('python contacts.py OUTPUT'.replace('OUTPUT', output_1))
            #os.system('python contacts.py OUTPUT'.replace('OUTPUT', output_2))

            contacts_1 = output_1[:-4] + '_contacts.csv'
            contacts_2 = output_2[:-4] + '_contacts.csv'

            os.system('python compare.py ' + contacts_1 + ' ' + contacts_2 + \
                      ' --flip1 --flip2')


#pdb_ids = []
#with open("../families_representatives.txt", "r", encoding="utf8") as f:
#    lines = f.readlines()
#    for line in lines:
#        if '>' not in line.strip():
#            pdb_ids.append(line.strip())
#pre_calculate(pdb_ids)
#exit(0)

pdb_ids = []
with open("../families_representatives.txt", "r", encoding="utf8") as f:
    lines = f.readlines()
    for _ in range(len(lines)):
        lines[_] = lines[_].replace("\n", "")
    family_path = lines[0].replace(">", "")
    print(family_path)
    for line in lines[1:]:
        if ">" in line:
            families_apd(pdb_ids, family_path)
            family_path = line[1:].replace(">", "")
            pdb_ids = []
        else:
            pdb_ids.append(line)
    families_apd(pdb_ids, family_path)

families_apd(pdb_ids, family_path)
