from Bio.Align import PairwiseAligner


def read_seq(seq_file_path: str) -> str:
    seq_file = open(seq_file_path, 'r')
    seq = seq_file.read().strip()
    seq_file.close()

    return seq


def read_chain(chain_file_path: str) -> str:
    chain_file = open(chain_file_path, 'r')
    residues = chain_file.readlines()
    chain_file.close()

    seq = ''

    resnum_ = None
    for residue in residues:
        resnum = int(residue.split()[-2])
        if resnum_ is not None:
            seq += '-' * (resnum - resnum_ - 1)
        seq += residue.split()[-1]
        resnum_ = resnum

    return seq


def align_seqs(seq_1: str, seq_2: str) -> tuple:
    aligner = PairwiseAligner(open_gap_score = -2, extend_gap_score = -1)
    alignment = aligner.align(seq_1, seq_2)[0]
    return alignment[0], alignment[1]


pdbid_1 = '6cu8'

for pdbid_2 in ['2n0a', '9kal', '7nck', '6peo', '9je7', '9je4', '8pk2',
                '6xyo', '9jdk', '9tpt']:
    seq_1 = read_seq('extracted_chains/' + pdbid_1 + '_seq.txt')
    seq_2 = read_seq('extracted_chains/' + pdbid_2 + '_seq.txt')

    chain_1 = read_chain('extracted_chains/' + pdbid_1 + '.txt')
    chain_2 = read_chain('extracted_chains/' + pdbid_2 + '.txt')

    #print('Sequences.')
    #print()
    #print(seq_1)
    #print()
    #print(seq_2)
    #print()

    #print('Chains.')
    #print()
    #print(chain_1)
    #print()
    #print(chain_2)
    #print()

    aligned_chain_1, aligned_chain_2 = align_seqs(seq_1, seq_2)

    print('Aligned sequences.')
    print()
    print(aligned_chain_1)
    print()
    print(aligned_chain_2)

    #print('Aligned sequences.')

    #seq_1_aligned, seq_2_aligned = align_seqs(seq_1, seq_2)

    #print()
    #print(seq_1_aligned)
    #print()
    #print(seq_2_aligned)
    #print()

    #seq_1_aligned, chain_1_aligned = align_seqs(seq_1, chain_1)

    #print('(1) Aligned sequence and chain.')
    #print()
    #print(seq_1_aligned)
    #print()
    #print(chain_1_aligned)
    #print()

    #seq_2_aligned, chain_2_aligned = align_seqs(seq_2, chain_2)

    #print('(2) Aligned sequence and chain.')
    #print()
    #print(seq_2_aligned)
    #print()
    #print(chain_2_aligned)
    #print()

    print('Aligned chains.')

    chain_1, chain_2 = align_seqs(chain_1, chain_2)

    print()
    print(chain_1)
    print()
    print(chain_2)
    print()

    #print('Resuling alignment of chains.')

    #
