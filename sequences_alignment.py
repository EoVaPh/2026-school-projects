from Bio import Align


def no_isolated_alignments(seq_1_aligned: str, seq_2_aligned: str) -> bool:
    '''Check whether the alignment represented by two strings of format AA--BB-
       AA---B- does not have isolated aligned residues like "-B-" in the
       example.'''

    assert len(seq_1_aligned) == len(seq_2_aligned)

    L = len(seq_1_aligned)

    for i in range(L):
        if seq_1_aligned[i] != '-':
            if i == 0 and seq_1_aligned[i+1] == '-':
                return False
            if i == L-1 and seq_1_aligned[i-1] == '-':
                return False
            if 0 < i < L-1 and seq_1_aligned[i-1] == '-' and \
                               seq_1_aligned[i+1] == '-':
                return False

        if seq_2_aligned[i] != '-':
            if i == 0 and seq_2_aligned[i+1] == '-':
                return False
            if i == L-1 and seq_2_aligned[i-1] == '-':
                return False
            if 0 < i < L-1 and seq_2_aligned[i-1] == '-' and \
                               seq_2_aligned[i+1] == '-':
                return False

    return True


def align_records(records_file_path_1: str, records_file_path_2: str) -> tuple:
    '''Find the best alignment of two given sequences among alignments which do
       not contain isolated aligned residues, and return the alignment in format
       of two aligned strings, that may contain "-" symbols for gaps.

    Parameters
    ----------
        records_file_path_1 : str
            Path to the first records file containing extended data about
            ATOMSEQ.
        records_file_path_2 : str
            Path to the second records file containing extended data about
            ATOMSEQ.

    Returns
    -------
        Two aligned sequences of the same length, where "-" symbols can be used
        for alignment with gaps.
    '''

    pars1 = []
    pars2 = []

    with open (records_file_path_1, "r", encoding = "utf8") as f:
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
    with open (records_file_path_2, "r", encoding = "utf8") as f:
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

    print('SEQUENCES TO ALIGN (INCLUDING UNSTRUCTURED)')
    print(seq1)
    print(seq2)

    aligner = Align.PairwiseAligner()
    aligner.match_score = 1.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -2.0

    alignments = aligner.align(seq1, seq2)

    for alignment in alignments:
        alignment_tokens = format(alignment, "phylip").split()
        seq_1_aligned, seq_2_aligned = alignment_tokens[2], alignment_tokens[3]
        print('ALIGNMENT')
        print(seq_1_aligned)
        print(seq_2_aligned)
        if no_isolated_alignments(seq_1_aligned, seq_2_aligned):
            print('ACCEPTED')
            best_alignment = alignment
            break

    print(f'Score: {best_alignment.score}.')

    alignment_tokens = format(best_alignment, "phylip").split()
    seq_1_aligned, seq_2_aligned = alignment_tokens[2], alignment_tokens[3]

    return seq_1_aligned, seq_2_aligned
