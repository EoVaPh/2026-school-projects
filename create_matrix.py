import numpy as np

amino_acids = "ARNDCQEGHILKMFPSTWYV"

match = 5
mismatch = -10

matrix = np.full(
    (len(amino_acids), len(amino_acids)),
    mismatch,
    dtype=int
)

np.fill_diagonal(matrix, match)

with open("matrix.txt", "w", encoding="utf-8") as f:

    f.write("   " + " ".join(amino_acids) + "\n")

    for aa, row in zip(amino_acids, matrix):
        f.write(
            f"{aa} "
            + " ".join(f"{value:3d}" for value in row)
            + "\n"
        )

print("matrix.txt created")