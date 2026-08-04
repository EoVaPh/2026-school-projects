from Bio import Align

pars1 = []
pars2 = []
with open ("records/1abc.txt", "r", encoding = "utf8") as f:
    preread1 = f.readlines()
    for i in preread1:
        i.replace("\n", "")
        i = i.split()
        i[0] = float(i[0])
        i[1] = float(i[1])
        i[2] = float(i[2])
        i[4] = int(i[4])
        pars1.append(i)
    print(pars1)
with open ("records/2abc.txt", "r", encoding = "utf8") as f:
    preread2 = f.readlines()
    for i in preread2:
        i.replace("\n", "")
        i = i.split()
        i[0] = float(i[0])
        i[1] = float(i[1])
        i[2] = float(i[2])
        i[4] = int(i[4])
        pars2.append(i)
    print(pars2)

seq1 = ""
seq2 = ""

res_id_now = 0
res_id_last = pars1[0][4]
for i in pars1:
    res_id_now = i[4]
    seq1 += ('-' * (res_id_now - res_id_last - 1))
    seq1 += (i[3])
    res_id_last = res_id_now

res_id_now = 0
res_id_last = pars2[0][4]
for i in pars2:
    res_id_now = i[4]
    seq2 += ("-" * (res_id_now - res_id_last - 1))
    seq2 += (i[3])
    res_id_last = res_id_now

print(seq1)
print(seq2)

aligner = Align.PairwiseAligner()
aligner.match_score = 1.0
aligner.mismatch_score = -1.0
aligner.open_gap_score = -2.0

alignments = aligner.align(seq1, seq2)
best_alignment = alignments[0]
print(f"Score: {best_alignment.score}\n")
print(best_alignment)