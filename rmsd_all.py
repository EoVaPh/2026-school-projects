from pathlib import Path
from itertools import combinations

import numpy as np
from scipy.spatial.transform import Rotation
from Bio.Align import PairwiseAligner


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


def get_positions(pdbid: str) -> list:
    '''Read atom coordinates from a chain file.'''

    chain_file = open('extracted_chains/' + f'{pdbid}.txt','r')
    residues = chain_file.readlines()
    chain_file.close()

    positions = []

    for residue in residues:
        residue_tokens = residue.strip().split()
        positions.append((
            float(residue_tokens[0]),
            float(residue_tokens[1]),
            float(residue_tokens[2])
        ))

    return positions


def align_seqs(seq_1: str, seq_2: str) -> tuple:
    '''Do global pairwise alignment of two amino acid sequences.'''

    aligner = PairwiseAligner(open_gap_score=-3,
                            extend_gap_score=-2,
                            mismatch_score=-1,
                            left_gap_score=-1,
                            right_gap_score=-1)
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
    '''Align a structural chain to its corresponding SEQRES
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
            chain_to_seqres[current_chain_position] = current_seqres_position

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
    '''Process one pair of structures:
    1. Align two SEQRES sequences.
    2. Map each SEQRES to the common alignment.
    3. Map each chain to its SEQRES.
    4. Project both chains onto the common SEQRES alignment.'''

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
    '''Remove leading and trailing positions containing only gaps.'''

    assert len(alignment_1) == len(alignment_2)

    N = len(alignment_1)

    for n in range(N):
        if (alignment_1[n] != '-' or alignment_2[n] != '-'):
            begin_n = n
            break

    for n in reversed(range(N)):
        if (alignment_1[n] != '-' or alignment_2[n] != '-'):
            end_n = n
            break

    return ''.join(alignment_1[begin_n:end_n + 1]),\
           ''.join(alignment_2[begin_n:end_n + 1])


def get_len_longest_shared_region(aligned_chain_1: str,
                                  aligned_chain_2: str) -> int:
    '''Find the length of the longest continuous region
    where the two chains contain the same amino acids.'''

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


def num_substitutions(seq_1: str, seq_2: str) -> int:
    '''Count amino acid substitutions between two equal-length sequences.'''

    assert len(seq_1) == len(seq_2)

    N = len(seq_1)

    count = 0

    for n in range(N):
        if seq_1[n] != seq_2[n]:
            count += 1

    return count


def get_shared_regions(w: int, aligned_chain_1: str, aligned_chain_2: str,
                               positions_1: list, positions_2: list) -> tuple:

    '''Find all identical continuous regions of length w.

    Windows containing gaps are excluded.
    Windows containing amino acid substitutions are excluded.'''

    assert len(aligned_chain_1) == len(aligned_chain_2)

    N = len(aligned_chain_1)

    i_1, i_2 = -1, -1

    aa_regions_1, pos_regions_1, aa_regions_2, pos_regions_2 = [], [], [], []

    for i in range(6, N - w + 1 - 6):
        if aligned_chain_1[i] != '-':
            i_1 += 1

        if aligned_chain_2[i] != '-':
            i_2 += 1

        region_sequence_1, region_positions_1 = [], []
        region_sequence_2, region_positions_2 = [], []

        for j in range(w):
            if aligned_chain_1[i + j] != '-' and aligned_chain_2[i + j] != '-':
                region_sequence_1 += aligned_chain_1[i + j]
                region_positions_1.append(positions_1[i_1 + j])
                region_sequence_2 += aligned_chain_2[i + j]
                region_positions_2.append(positions_2[i_2 + j])
            else:
                break

        if len(region_sequence_1) == w:
            if num_substitutions(region_sequence_1, region_sequence_2) <= 0:
                aa_regions_1.append(region_sequence_1)
                pos_regions_1.append(region_positions_1)
                aa_regions_2.append(region_sequence_2)
                pos_regions_2.append(region_positions_2)

    return aa_regions_1, pos_regions_1, aa_regions_2, pos_regions_2


def calc_rmsd(positions_1: list, positions_2: list) -> float:
    '''Calculate RMSD between two lists of 3D points.'''

    assert len(positions_1) == len(positions_2)

    positions_1_array = np.array(positions_1)
    positions_2_array = np.array(positions_2)

    centroid_1 = np.mean(positions_1_array, axis=0)
    centroid_2 = np.mean(positions_2_array, axis=0)

    positions_1_centered = (positions_1_array - centroid_1)
    positions_2_centered = (positions_2_array - centroid_2)

    rotation, rmsd_value = Rotation.align_vectors(
        positions_1_centered, positions_2_centered
    )

    return rmsd_value


input_folder = Path('extracted_chains')
output_file = Path('rmsds.txt')

pdb_ids = sorted(f.name[:-8] for f in input_folder.glob('*_seq.txt') if f.is_file())
print(f'Found {len(pdb_ids)} structures.')

pairs = list(combinations(pdb_ids, 2))
print(f'Total unique pairs: {len(pairs)}')

with open(output_file, 'w') as output:

    for pair_number, (pdbid_1, pdbid_2) in enumerate(pairs, start=1):

        print(
            f'[{pair_number}/{len(pairs)}] '
            f'{pdbid_1} / {pdbid_2}'
        )

        try:
            results = process_pair(pdbid_1, pdbid_2)

            aligned_chain_1, aligned_chain_2 = strip_alignment(
                results['chain_alignment'][0], results['chain_alignment'][1]
            )

            len_longest_shared_region = get_len_longest_shared_region(
                aligned_chain_1, aligned_chain_2
            )

            if len_longest_shared_region >= 7:

                positions_1 = get_positions(pdbid_1)
                positions_2 = get_positions(pdbid_2)

                aa_regions_1, pos_regions_1, aa_regions_2, pos_regions_2 = \
                    get_shared_regions(
                        6,
                        aligned_chain_1, aligned_chain_2,
                        positions_1, positions_2
                    )

                rmsd_length_6 = []

                for r in range(len(pos_regions_1)):
                    rmsd_length_6.append(
                        calc_rmsd(pos_regions_1[r], pos_regions_2[r])
                    )

                if rmsd_length_6:

                    output.write(
                        f'>{pdbid_1}_{pdbid_2}\n'
                    )

                    output.write(
                        ' '.join(
                            f'{value:.3f}'
                            for value in rmsd_length_6
                        )
                        + '\n'
                    )

        except Exception as e:
            print(f'ERROR: {pdbid_1} / {pdbid_2}: {e}')


print()
print('Finished.')
print(f'Results saved to: {output_file}')