from pathlib import Path

amino_acids = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V"
}

def get_amino_acid(code):
    return amino_acids.get(code, "X")


folder = Path(r"C:\Users\User\Documents\bioinf\smtb\all_amyloid_structures\amyloid_structures")

for structure_file in folder.glob("*.pdb"):
    name = structure_file.name
    print(f"Обрабатываю: {structure_file.name}")

    with open(structure_file, "r", encoding="utf-8") as file:
        structure_lines = file.readlines()

    records = []
    seq = ''
    chain = ''
    seq_res = ''

    for line in structure_lines:
        if 'SEQRES' == line[:6]:
            elem = line.split()
            if chain == '':
                chain = elem[2]
            elif chain != elem[2]:
                continue
            for aminoacid in elem[4:]:
                seq_res += get_amino_acid(aminoacid)

        elif ('ATOM  ' == line[:6] or 'HETATM' == line[:6]) and ' CA ' == line[12:16] and line[16] in (' ', 'A'):
            one_record = (float(line[30:38]), float(line[38:46]), float(line[46:54]), get_amino_acid(line[17:20]), int(line[22:26]))
        #    one_record = [line[12:16], line[16:17], amino_acids[line[17:20]], line[21:22], int(line[22:26]), line[26:27], (float(line[30:38]), float(line[38:46]), float(line[46:54]))]
        #    one_record = [line[:6], line[6:11], line[12:16], line[16:17], line[17:20], line[21:22], line[22:26], line[26:27], line[30:38], line[38:46], line[46:54], line[54:60], line[60:66], line[76:78], line[78:80]]
            records.append(one_record)
            seq += get_amino_acid(line[17:20])
            chain = line[21:22]
        elif 'TER   ' == line[:6]:
            break


    with open(folder.parent / f"{name}_sequences.fa", "w") as file:
        print(">SEQRES\n", seq_res, file=file)
        print(">ATOM sequence\n", seq, file=file)
    with open(folder.parent / f"{name}_records.txt", "w") as file:
        for record in records:
            for item in record:
                print(item, end='\t', file=file)
            print('', file=file)