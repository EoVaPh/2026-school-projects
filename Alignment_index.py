input_file = "alignment.txt"
output_file = "pairs.txt"

with open(r'C:\Users\User\Documents\bioinf\smtb\alignment.txt', "r", encoding="utf-8") as file:
    alignment = file.readlines()

seq1 = alignment[0].strip()
seq2 = alignment[1].strip()

if len(seq1) != len(seq2):
    raise ValueError("Aligned sequences must have the same length")

pairs = []

index1 = 0
index2 = 0

for aa1, aa2 in zip(seq1, seq2):

    if aa1 != "-" and aa2 != "-":
        pairs.append((index1, index2))

    if aa1 != "-":
        index1 += 1

    if aa2 != "-":
        index2 += 1


with open('alignment_index.txt', "w", encoding="utf-8") as file:
    for index1, index2 in pairs:
        print(index1, index2, file=file)
