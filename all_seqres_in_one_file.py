from pathlib import Path


input_folder = Path(
    r"C:\Users\User\Documents\bioinf\smtb\sequences"
)

output_folder = Path(
    r"C:\Users\User\Documents\bioinf\smtb"
)

output_folder.mkdir(exist_ok=True)

output_file = output_folder / "seqres_all.txt"


with open(output_file, "w", encoding="utf8") as output:

    for file in input_folder.glob("*_sequences.fa"):

        with open(file, "r", encoding="utf8") as input_file:
            lines = input_file.readlines()

        if len(lines) < 2:
            continue

        sequence = lines[1].strip()

        file_name = file.stem[:4]

        print(f">{file_name}", file=output)
        print(sequence, file=output)
