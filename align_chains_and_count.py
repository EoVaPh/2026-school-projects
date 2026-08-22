from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation

from Bio.Align import PairwiseAligner

from aiupred import AIUPred

import math
from typing import List, Tuple, Sequence

from matplotlib import pyplot as plot
import matplotlib

matplotlib.rcParams['figure.dpi'] = 300
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rc('font', family='STIXGeneral')
matplotlib.rc('font', weight='ultralight')


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

    aligner = PairwiseAligner(open_gap_score = -3,
                              extend_gap_score = -2,
                              mismatch_score = -1,
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


def get_positions(pdbid: str) -> list:
    chain_file = open('extracted_chains/' + pdbid + '.txt', 'r')
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


def num_substitutions(seq_1: str, seq_2: str) -> int:
    assert len(seq_1) == len(seq_2)

    N = len(seq_1)

    count = 0

    for n in range(N):
        if seq_1[n] != seq_2[n]:
            count += 1

    return count


def get_shared_regions(w: int, aligned_chain_1: str, aligned_chain_2: str,
                               positions_1: list, positions_2: list) -> tuple:
    '''Find all shared continuous regions of the same length and return their
       sequences and coordinates.'''

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
            if num_substitutions(region_sequence_1, region_sequence_2) <= 1:
                aa_regions_1.append(region_sequence_1)
                pos_regions_1.append(region_positions_1)
                aa_regions_2.append(region_sequence_2)
                pos_regions_2.append(region_positions_2)

    return aa_regions_1, pos_regions_1, aa_regions_2, pos_regions_2


def calc_squared_distance(pos_1: tuple, pos_2: tuple) -> float:
    return (pos_1[0] - pos_2[0])**2 + \
           (pos_1[1] - pos_2[1])**2 + \
           (pos_1[2] - pos_2[2])**2


def calc_rmsd(positions_1: list, positions_2: list) -> float:
    '''Calculate RMSD of two lists of points.'''

    assert len(positions_1) == len(positions_2)

    positions_1_array = np.array(positions_1)
    positions_2_array = np.array(positions_2)

    centroid_1 = np.mean(positions_1_array, axis=0)
    centroid_2 = np.mean(positions_2_array, axis=0)

    positions_1_centered = positions_1_array - centroid_1
    positions_2_centered = positions_2_array - centroid_2

    rotation, rmsd_value = Rotation.align_vectors(
        positions_1_centered, positions_2_centered
    )

    return rmsd_value


def calc_lddt(
    reference: List[Tuple[float, float, float]],
    model: List[Tuple[float, float, float]],
    inclusion_radius: float = 15.0,
    thresholds: Sequence[float] = (0.5, 1.0, 2.0, 4.0)
) -> float:
    """
    Calculate the Local Distance Difference Test (LDDT) score between two
    structures represented as lists of corresponding 3D coordinates.

    Parameters
    ----------
    reference : list of tuples
        Reference coordinates.
    model : list of tuples
        Model coordinates, same order and length as reference.
    inclusion_radius : float
        Include reference pairs whose distance is <= this value.
    thresholds : sequence of floats
        Distance-difference thresholds used by LDDT.

    Returns
    -------
    float
        The LDDT score.
    """
    if len(reference) != len(model):
        raise ValueError("reference and model must have the same length")
    if len(reference) < 2:
        raise ValueError("at least two points are required")
    if not thresholds:
        raise ValueError("thresholds must not be empty")

    def distance(a: Tuple[float, float, float],
                 b: Tuple[float, float, float]) -> float:
        return math.sqrt(
            (a[0] - b[0]) ** 2
            + (a[1] - b[1]) ** 2
            + (a[2] - b[2]) ** 2
        )

    n = len(reference)
    total_score = 0.0
    considered_points = 0

    for i in range(n):
        ref_i = reference[i]
        model_i = model[i]
        neighbor_count = 0
        preserved_count = 0

        for j in range(n):
            if i == j:
                continue

            d_ref = distance(ref_i, reference[j])
            if d_ref > inclusion_radius:
                continue

            d_model = distance(model_i, model[j])
            diff = abs(d_model - d_ref)

            for threshold in thresholds:
                if diff < threshold:
                    preserved_count += 1

            neighbor_count += 1

        # Ignore points with no neighbours inside the inclusion radius.
        if neighbor_count == 0:
            continue

        point_score = preserved_count / (neighbor_count * len(thresholds))
        total_score += point_score
        considered_points += 1

    if considered_points == 0:
        return 0.0

    return total_score / considered_points


predictor = AIUPred()

def calculate_seqres_idr(pdb_id: str, folder: str) -> list:
    '''Find a SEQRES file by PDB ID, calculate AIUPred disorder
    propensities, and return them as a list.'''

    folder = Path(folder)
    file_path = folder / f'{pdb_id}_seq.txt'

    if not file_path.exists():
        raise FileNotFoundError(f'File not found: {file_path}')

    with open(file_path, 'r') as file:
        sequence = file.read().strip()

    disorder_propensities = predictor.predict_disorder(sequence)

    return disorder_propensities.tolist()


def select_atomseq_idr(seqres_aligned: str, chain_aligned: str, aiupred_values: list) -> list:
    '''Having pre-calculated IDR propensities for a SEQRES, transfer them to the
    chain in accordance with SEQRES-chain alignment.'''

    assert len(seqres_aligned) == len(chain_aligned)

    IDR = []

    seqres_index = 0

    for i in range(len(seqres_aligned)):
        if seqres_aligned[i] != "-":
            p = aiupred_values[seqres_index]
            seqres_index += 1
        if chain_aligned[i] != "-":
            IDR.append(p)

    return IDR


def find_mismatches(aligned_chain_1: str, aligned_chain_2: str) -> str:
    assert len(aligned_chain_1) == len(aligned_chain_2)

    N = len(aligned_chain_1)

    mismatches = ''

    for n in range(N):
        if aligned_chain_1[n] != '-' and aligned_chain_2[n] != '-' and \
           aligned_chain_1[n] != aligned_chain_2[n]:
                mismatches += '*'
        else:
            mismatches += ' '

    return mismatches


pdb_ids = [f.name[:4] for f in Path('extracted_chains').rglob('*') \
           if f.is_file() and '.txt' in f.name]

verbose = False

rmsd_length_6, lddt_length_6, idr_length_6 = [], [], []

pdb_id = '8pkg'

for pdbid_1 in [pdb_id]:#pdb_ids:
    for pdbid_2 in pdb_ids:
        results = process_pair(pdbid_1, pdbid_2)

        aligned_chain_1, aligned_chain_2 = strip_alignment(
            results['chain_alignment'][0], results['chain_alignment'][1]
        )

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
            print(find_mismatches(aligned_chain_1, aligned_chain_2))
            print()

            positions_1 = get_positions(pdbid_1)
            positions_2 = get_positions(pdbid_2)

            aa_regions_1, pos_regions_1, aa_regions_2, pos_regions_2 = \
                get_shared_regions(
                    6,
                    aligned_chain_1, aligned_chain_2,
                    positions_1, positions_2
                )

            seqres_idr_1 = calculate_seqres_idr(pdbid_1, 'extracted_chains')
            seqres_idr_2 = calculate_seqres_idr(pdbid_2, 'extracted_chains')

            chain_idr_1 = select_atomseq_idr(
                results['chain_1'][0],
                results['chain_1'][1],
                seqres_idr_1
            )

            chain_idr_2 = select_atomseq_idr(
                results['chain_2'][0],
                results['chain_2'][1],
                seqres_idr_2
            )

            idr_by_position_1 = dict(zip(positions_1, chain_idr_1))
            idr_by_position_2 = dict(zip(positions_2, chain_idr_2))

            for r in range(len(pos_regions_1)):
                rmsd_length_6.append(
                    calc_rmsd(pos_regions_1[r], pos_regions_2[r])
                )

                lddt_length_6.append(
                    1 - calc_lddt(pos_regions_1[r], pos_regions_2[r])
                )

                window_pair_idr = []

                for pos_1, pos_2 in zip(pos_regions_1[r], pos_regions_2[r]):

                    idr_1 = idr_by_position_1[pos_1]
                    idr_2 = idr_by_position_2[pos_2]
                    mean_idr_residue = (idr_1 + idr_2) / 2
                    window_pair_idr.append(mean_idr_residue)

                mean_idr_window = np.mean(window_pair_idr)
                idr_length_6.append(mean_idr_window)

            #num_residues_1 = len(aligned_chain_1) - aligned_chain_1.count('-')
            #num_residues_2 = len(aligned_chain_2) - aligned_chain_2.count('-')


plot.scatter(idr_length_6, lddt_length_6, color='#a53860', alpha=0.1)
plot.xticks(fontsize=12)
plot.yticks(fontsize=12)
plot.xlabel('IDR', fontsize=16)
#plot.ylabel('RMSD', fontsize=16)
plot.ylabel('1 – LDDT', fontsize=16)
plot.title(pdb_id + ' family screened using window of length 6', fontsize=16)
plot.tight_layout()
plot.savefig('idr_lddt.png')
