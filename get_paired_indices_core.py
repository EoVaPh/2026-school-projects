def get_paired_indices_core(seq_1: str, seq_2: str,
                            core_begin: int, core_end: int) -> list:
    N = len(seq_1)

    assert len(seq_2) == N

    i, i_1, i_2 = 0, -1, -1
    paired_indices_core = []

    for i in range(N):
        if seq_1[i] != '-':
            i_1 += 1

        if seq_2[i] != '-':
            i_2 += 1

        if seq_1[i] != '-' and seq_2[i] != '-' and core_begin <= i <= core_end:
            paired_indices_core.append((i_1, i_2))

    return paired_indices_core


#print(get_paired_indices_core('---AAA-AA-A-',
#                              '--AAAAAAAAAA',
#                              4, 10))
