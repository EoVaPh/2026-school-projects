from Bio.Align import PairwiseAligner
from pathlib import Path


def read_seq(seq_file_path: str) -> str:
    '''Read a SEQRES sequence from a text file.'''

    seq_file = open(seq_file_path, 'r')
    seq = seq_file.read().strip()
    seq_file.close()

    return seq


def read_chain(chain_file_path: str) -> str:
    '''Read a structural chain from a text file.'''

    chain_file = open(chain_file_path, 'r')
    residues = chain_file.readlines()
    chain_file.close()

    seq = ''

    for residue in residues:
        resnum = int(residue.split()[-2])
        seq += residue.split()[-1]

    return seq


def align_seqs(seq_1: str, seq_2: str) -> tuple:
    '''Do global pairwise alignment of two amino acid sequences.'''

    aligner = PairwiseAligner(open_gap_score = -2,
                              extend_gap_score = -1,
                              left_gap_score = -1,
                              right_gap_score = -1)
    alignment = aligner.align(seq_1, seq_2)[0]
    return alignment[0], alignment[1]


def mapping_seqres_to_common_alignment(aligned_seq: str) -> dict:
    '''Map SEQRES positions to positions in the common alignment.'''

    seqres_to_common_alignment = {}
    seqres_position = 0

    for alignment_position, aa in enumerate(aligned_seq):
        if aa != '-':
            seqres_to_common_alignment[seqres_position] = alignment_position
            seqres_position += 1

    return seqres_to_common_alignment


def mapping_chain_to_seqres(seqres: str, chain: str) -> tuple:
    '''Align a structural chain to its corresponding SEQRES sequence
    and map chain positions to SEQRES positions.'''

    aligned_seqres, aligned_chain = align_seqs(seqres, chain)

    chain_to_seqres = {}
    seqres_position = 0
    chain_position = 0

    for aa_seqres, aa_chain in zip(aligned_seqres, aligned_chain):
        current_seqres_position = None
        current_chain_position = None

        if aa_seqres != '-':
            current_seqres_position = seqres_position
            seqres_position += 1

        if aa_chain != '-':
            current_chain_position = chain_position
            chain_position += 1

        if (current_seqres_position is not None and current_chain_position is not None):
            chain_to_seqres[current_chain_position] = (current_seqres_position)

    return aligned_seqres, aligned_chain, chain_to_seqres


def project_chain_to_common_alignment(chain: str, chain_to_seqres: dict, seqres_to_common_alignment: dict, alignment_length: int) -> list:
    '''Project a structural chain onto the common SEQRES alignment.'''

    result = ['-'] * alignment_length

    for chain_position, seqres_position in chain_to_seqres.items():
        if seqres_position not in seqres_to_common_alignment:
            continue

        alignment_position = seqres_to_common_alignment[seqres_position]

        result[alignment_position] = chain[chain_position]

    return result


def process_pair(pdbid_1: str, pdbid_2: str) -> dict:
    '''Process a pair of structures:
    1. Align two SEQRES sequences
    2. Map each SEQRES to the common alignment
    3. Map each chain to it's SEQRES
    4. Project both chains onto the common SEQRES alignment'''

    seq_1 = read_seq('extracted_chains/' + pdbid_1 + '_seq.txt')
    seq_2 = read_seq('extracted_chains/' + pdbid_2 + '_seq.txt')

    chain_1 = read_chain('extracted_chains/' + pdbid_1 + '.txt')
    chain_2 = read_chain('extracted_chains/' + pdbid_2 + '.txt')

    aligned_seqres_1, aligned_seqres_2 = align_seqs(seq_1, seq_2)

    seqres_1_to_common = mapping_seqres_to_common_alignment(aligned_seqres_1)
    seqres_2_to_common = mapping_seqres_to_common_alignment(aligned_seqres_2)
    common_alignment_length = len(aligned_seqres_1)

    chain_1_seqres_alignment, chain_1_chain_alignment, chain_1_to_seqres = mapping_chain_to_seqres(seq_1, chain_1)

    chain_2_seqres_alignment, chain_2_chain_alignment, chain_2_to_seqres = mapping_chain_to_seqres(seq_2, chain_2)

    projected_chain_1 = project_chain_to_common_alignment(chain_1, chain_1_to_seqres, seqres_1_to_common, common_alignment_length)
    projected_chain_2 = project_chain_to_common_alignment(chain_2, chain_2_to_seqres, seqres_2_to_common, common_alignment_length)

    return {
        'seqres_alignment': (
            aligned_seqres_1,
            aligned_seqres_2
        ),

        'chain_1': (
            chain_1_seqres_alignment,
            chain_1_chain_alignment
        ),

        'chain_2': (
            chain_2_seqres_alignment,
            chain_2_chain_alignment
        ),

        'chain_alignment': (
            projected_chain_1,
            projected_chain_2
        ),

        'chain_1_to_seqres': chain_1_to_seqres,
        'chain_2_to_seqres': chain_2_to_seqres
    }


def strip_alignment(alignment_1, alignment_2) -> tuple:
    assert len(alignment_1) == len(alignment_2)

    N = len(alignment_1)

    for n in range(0, N):
        if alignment_1[n] != '-' or alignment_2[n] != '-':
            begin_n = n
            break

    for n in reversed(range(0, N)):
        if alignment_1[n] != '-' or alignment_2[n] != '-':
            end_n = n
            break

    return ''.join(alignment_1[begin_n:end_n+1]),\
           ''.join(alignment_2[begin_n:end_n+1])


def get_len_longest_shared_region(aligned_chain_1: str,
                                  aligned_chain_2: str) -> int:
    assert len(aligned_chain_1) == len(aligned_chain_2)

    len_shared_region, len_longest_shared_region = 0, 0

    N = len(aligned_chain_1)

    for n in range(N):
        if aligned_chain_1[n] != '-':
            if aligned_chain_1[n] == aligned_chain_2[n]:
                len_shared_region += 1
                len_longest_shared_region = max(len_longest_shared_region,
                                                len_shared_region)
            else:
                len_shared_region = 0

    return len_longest_shared_region


pdbid_1 = '9k23'

pdb_ids = [f.name[:4] for f in Path('extracted_chains').rglob('*') if f.is_file() and '.txt' in f.name]

verbose = False

for pdbid_2 in pdb_ids:
    results = process_pair(pdbid_1, pdbid_2)

    aligned_chain_1, aligned_chain_2 = strip_alignment(
        results['chain_alignment'][0], results['chain_alignment'][1]
    )

    #num_residues_1 = len(aligned_chain_1) - aligned_chain_1.count('-')
    #num_residues_2 = len(aligned_chain_2) - aligned_chain_2.count('-')

    #print(num_residues_1, num_residues_2)

    len_longest_shared_region = get_len_longest_shared_region(
        aligned_chain_1, aligned_chain_2
    )

    if len_longest_shared_region >= 7:
        print(f'{pdbid_1} / {pdbid_2}')

        if verbose:
            print("SEQRES alignment")
            print()
            print(results['seqres_alignment'][0])
            print(results['seqres_alignment'][1])
            print()

            print("CHAIN 1 -> SEQRES 1")
            print()
            print(results['chain_1'][0])
            print(results['chain_1'][1])
            print()

            print("CHAIN 2 -> SEQRES 2")
            print()
            print(results['chain_2'][0])
            print(results['chain_2'][1])
            print()

            print("FINAL CHAIN ALIGNMENT")

        print()
        print(aligned_chain_1)
        print(aligned_chain_2)
        print()
