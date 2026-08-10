from pathlib import Path
from Bio import SeqIO


fasta_file = Path("seqres_all.txt")
families_file = Path("families.txt")

output_folder = Path("families_seqres")
output_folder.mkdir(exist_ok=True)


families = {}

current_family = None

with open(families_file, "r", encoding="utf-8") as file:

    for line in file:
        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):
            current_family = line[1:].strip()
            families[current_family] = []

        else:
            families[current_family].append(line)

records = list(SeqIO.parse(fasta_file, "fasta"))


for number, (family_name, structures) in enumerate(
    families.items(), start=1
):

    family_code = f"family_{number:03d}"

    structures = set(structures)

    selected_records = []

    for record in records:

        structure_name = record.id[:4]

        if structure_name in structures:
            selected_records.append(record)


    output_file = output_folder / f"{family_code}.fasta"

    SeqIO.write(
        selected_records,
        output_file,
        "fasta"
    )

    print(
        f"{family_code}: "
        f"{len(selected_records)} sequences"
    )