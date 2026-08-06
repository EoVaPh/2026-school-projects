from Bio import Align
import math

bypair = []
def aligners (path1, path2):

    # collecting sequences
    with open (path1, "r", encoding="utf8") as f:
        f.readline()
        seq1 = f.readline().strip()
        if not seq1:
            seq1 = f.readline().strip()
    with open (path2, "r", encoding="utf8") as f:
        f.readline()
        seq2 = f.readline().strip()
        if not seq2:
            seq2 = f.readline().strip()
    # rename pathways to pdbid of amyloid and make a list with ids
    path1 = str(path1)
    path1 = path1[30:34]
    path2 = str(path2)
    path2 = path2[30:34]
    pair = []
    pair.append(path1)
    pair.append(path2)

    #make an alignment
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 1.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -2.0
    alignments = aligner.align(seq1, seq2)
    score_AA = aligner.score(seq1, seq1)
    score_BB = aligner.score(seq2, seq2)
    #normalizing a score
    if alignments:
        alignment = alignments[0]
        score_raw = alignment.score
        normalized_score = score_raw/(math.sqrt(score_AA*score_BB))
        pair.append(normalized_score)
    else:
        pair.append(0)
    bypair.append(pair)
    return(bypair)