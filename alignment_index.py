def get_aligned_indices(seq_1_aligned: str, seq_2_aligned: str):
    seq1 = seq_1_aligned.strip()
    seq2 = seq_2_aligned.strip()

    if len(seq1) != len(seq2):
        raise ValueError('Aligned sequences must have the same length.')

    pairs = []

    index1 = 0
    index2 = 0

    for aa1, aa2 in zip(seq1, seq2):
        if aa1 != '-' and aa2 != '-':
            pairs.append((index1, index2))

        if aa1 != '-':
            index1 += 1

        if aa2 != '-':
            index2 += 1

    return pairs
