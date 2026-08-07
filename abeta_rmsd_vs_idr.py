import numpy as np
from matplotlib import pyplot as plot
import matplotlib
import os
from itertools import combinations
from contextlib import redirect_stdout
from io import StringIO


from pathlib import Path

from Bio import Align

from sequences_alignment import align_records
from sequences_alignment import no_isolated_alignments
from alignment_index import get_aligned_indices

import plotly.tools as tools

matplotlib.rcParams['figure.dpi'] = 300
matplotlib.rcParams['mathtext.fontset'] = 'stix'
#matplotlib.rc('font', family='STIXGeneral')
#matplotlib.rc('font', weight='ultralight')


def get_plane_basis(points_1: list, points_2: list, paired_indices: list):
    if paired_indices is None:
        points = points_1 + points_2
    else:
        points = []
        for pair in paired_indices:
            i, j = pair[0], pair[1]
            points.append(points_1[i])
            points.append(points_2[j])

    pts = np.array(points, dtype=float)
    mean = np.mean(pts, axis=0)
    centered = pts - mean

    cov = np.cov(centered, rowvar=False)

    eigvals, eigvecs = np.linalg.eigh(cov)

    plane_basis = eigvecs[:, 1:]

    return plane_basis


def draw_projection(points: list, res_ids: list, seq: str, color: str,
                    zorder: int, plane_basis):
    '''Draw projection of the points on a plane with defined basis.'''

    pts = np.array(points, dtype=float)
    mean = np.mean(pts, axis=0)
    centered = pts - mean

    coords_2d = centered @ plane_basis

    projections = [tuple(coord) for coord in coords_2d]

    xs, ys = [coord[0] for coord in projections], [coord[1] for coord in projections]

    plot.scatter(xs, ys, color=color, s=256, zorder=zorder)

    for i in range(len(projections)):
        plot.text(xs[i] - 0.25, ys[i] - 0.35, seq[i], color='white',
                  zorder=zorder, fontweight='bold')

    for i in range(len(points) - 1):
        if res_ids[i+1] - res_ids[i] == 1:
            plot.plot([xs[i], xs[i+1]], [ys[i], ys[i+1]], color=color,
                      linewidth=4, zorder=zorder-1, alpha=0.5)

    return plane_basis, xs, ys


def read_records(records_file_path: str):
    '''Read points, one-letter sequence, and residue numbers from a records
       file.'''

    records_file = open(records_file_path, 'r')
    records = records_file.readlines()
    records_file.close()

    points, seq, res_ids = [], '', []

    for record in records:
        tokens = record.split()
        points.append((float(tokens[0]), float(tokens[1]), float(tokens[2])))
        seq += tokens[3]
        res_ids.append(int(tokens[4]))

    return points, seq, res_ids


def read_seqres(file_path: str, structure_name: str) -> str:
    '''Read SEQRES from a file containing SEQRES of the structures, in
       FASTA format.

       >PDB ID
       SEQRES'''

    with open(file_path, "r", encoding="utf8") as file:
        lines = file.readlines()
    for i in range(len(lines)):
        if lines[i].strip() == ">" + structure_name:
            return lines[i + 1].strip()

    raise ValueError(f"SEQRES for {structure_name} not found")


def read_aiupred(file_path: str, structure_name: str) -> list:
    '''Read IDR propensities of a SEQRES from a FASTA-like file in format

    >PDB ID
    p_1 p_2 p_3'''

    with open(file_path, "r", encoding="utf8") as file:
        lines = file.readlines()

    for i in range(len(lines)):
        if lines[i].strip() == ">" + structure_name:
            values = lines[i + 1].split()
            return [float(value) for value in values]

    raise ValueError(f"AIUPred values for {structure_name} not found")


def align_sequences(seq_1, seq_2):

    aligner = Align.PairwiseAligner()

    aligner.match_score = 1.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -2.0

    alignments = aligner.align(seq_1, seq_2)

    for alignment in alignments:

        tokens = format(alignment, "phylip").split()

        seq_1_aligned = tokens[2]
        seq_2_aligned = tokens[3]

        if no_isolated_alignments(seq_1_aligned, seq_2_aligned):
            return (seq_1_aligned, seq_2_aligned)

    raise ValueError("No suitable alignment found")


def select_atomseq_idr(seqres: str, atom_seq: str, aiupred_values: list) -> list:
    '''Having pre-calculated IDR propensities for a SEQRES, transit them ATOMSEQ
       in accordance with SEQRES-ATOMSEQ alignment.'''

    seqres_aligned, atom_aligned = align_sequences(seqres, atom_seq)

    print("\nSEQRES / ATOM alignment:")
    print("SEQRES:", seqres_aligned)
    print("ATOM:  ", atom_aligned)

    IDR = []

    seqres_index = 0

    for i in range(len(seqres_aligned)):
        if seqres_aligned[i] != "-":
            p = aiupred_values[seqres_index]
            seqres_index += 1
        if atom_aligned[i] != "-":
            IDR.append(p)

    return IDR


def best_rmsd_transform(list1: list, list2: list, pairs: list) -> list:
    """
    Returns a new list of points that is the second list transformed by the
    optimal rigid transformation (rotation + translation) to minimise the
    Root Mean Square Deviation (RMSD) for the paired points.

    Parameters
    ----------
    list1 : list of tuple of (float, float, float)
        Reference points.
    list2 : list of tuple of (float, float, float)
        Points to be transformed.
    pairs : list of tuple of (int, int)
        Correspondences as (index_in_list1, index_in_list2).

    Returns
    -------
    transformed : list of tuple of (float, float, float)
        All points from list2 after applying the optimal rotation and
        translation.
    """
    # Convert to numpy arrays (Nx3)
    ref = np.array(list1, dtype=float)
    mob = np.array(list2, dtype=float)

    if len(pairs) == 0:
        raise ValueError("At least one pair is required to define a transformation.")

    # Extract paired points using the given indices
    idx1 = [i for i, j in pairs]
    idx2 = [j for i, j in pairs]
    P = ref[idx1]   # reference paired points
    Q = mob[idx2]   # mobile paired points

    # Centroids of the paired subsets
    c_ref = P.mean(axis=0)
    c_mob = Q.mean(axis=0)

    # Centre the paired points
    P_centered = P - c_ref
    Q_centered = Q - c_mob

    # Covariance matrix: H = Q_centered.T @ P_centered
    H = Q_centered.T @ P_centered

    # Singular Value Decomposition
    U, S, Vt = np.linalg.svd(H)

    # Optimal rotation (as a matrix that multiplies a row vector on the right)
    R = U @ Vt

    # Ensure a proper rotation (det = +1), i.e. no reflection.
    # If your use case allows reflections, simply remove this block.
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = U @ Vt

    # Apply the transformation to the whole second list:
    #   transformed = (mob - c_mob) @ R + c_ref
    transformed = (mob - c_mob) @ R + c_ref

    # Convert back to list of tuples
    return [tuple(pt) for pt in transformed]


def get_window_rmsd(seq_1_aligned, seq_2_aligned, points_1, points_2, window_size):


    paired_indices = get_aligned_indices(seq_1_aligned, seq_2_aligned)

    results = {}

    for start in range(len(paired_indices) - window_size + 1):

        end = start + window_size

        frame_name = f"{start} - {end - 1}"

        frame_pairs = paired_indices[start:end]

        first_pair = frame_pairs[0]
        last_pair = frame_pairs[-1]

        index1 = 0
        index2 = 0

        first_position = None
        last_position = None

        for position, (aa1, aa2) in enumerate(zip(seq_1_aligned, seq_2_aligned)):

            if aa1 != '-' and aa2 != '-':
                current_pair = (index1, index2)

                if current_pair == first_pair:
                    first_position = position

                if current_pair == last_pair:
                    last_position = position

            if aa1 != '-':
                index1 += 1

            if aa2 != '-':
                index2 += 1

        window_1 = seq_1_aligned[first_position:last_position + 1]
        window_2 = seq_2_aligned[first_position:last_position + 1]

        if '-' in window_1 or '-' in window_2:
            results[frame_name] = None
            continue

        frame_points_1 = [points_1[i] for i, j in frame_pairs]
        frame_points_2 = [points_2[j] for i, j in frame_pairs]

        local_pairs = [(i, i) for i in range(window_size)]
        transformed_points_2 = best_rmsd_transform(frame_points_1, frame_points_2, local_pairs)

        squared_distances = []

        for p1, p2 in zip(frame_points_1, transformed_points_2):
            distance_squared = np.sum((np.array(p1) - np.array(p2)) ** 2)
            squared_distances.append(distance_squared)

        rmsd = np.sqrt(np.mean(squared_distances))

        results[frame_name] = rmsd

    return results


def plot_rmsd_vs_mean_idr(results: dict, IDR_1: list, IDR_2: list):


    mean_idr_values = []
    rmsd_values = []

    for frame_name, rmsd in results.items():

        if rmsd is None:
            continue

        start, end = map(int, frame_name.split(" - "))

        idr_values = []

        for pair_index in range(start, end + 1):

            idr_1 = IDR_1[pair_index]
            idr_2 = IDR_2[pair_index]

            mean_pair_idr = (idr_1 + idr_2) / 2
            idr_values.append(mean_pair_idr)

        mean_idr = np.mean(idr_values)
        mean_idr_values.append(mean_idr)
        rmsd_values.append(rmsd)

    return mean_idr_values, rmsd_values



def connect_projected_points(xs_1: list, ys_1: list, xs_2: list, ys_2: list,
                             paired_indices: list):
    for pair in paired_indices:
        i, j = pair[0], pair[1]
        plot.plot([xs_1[i], xs_2[j]], [ys_1[i], ys_2[j]], linewidth=1,
                  color='black', zorder = -1, linestyle='--')


def show_points_matplotlib(list1, color1, name1, IDR_1, IDR_1_value,
                           list2, color2, name2, IDR_2, IDR_2_value):
    fig = plot.figure()
    ax = fig.add_subplot(111, projection='3d')

    if list1:
        x1, y1, z1 = zip(*list1)
        ax.plot(x1, y1, z1, color=color1, linewidth=1, label=name1)

        idr_indices_1 = [i for i, value in enumerate(IDR_1) if value > IDR_1_value]
        print(f"\n{name1}: atoms with IDR > {IDR_1_value}:")
        print(idr_indices_1)

        idr_points_1 = [list1[i] for i in idr_indices_1 if i < len(list1)]
        idr_points_1 = idr_points_1[:]

        if idr_points_1:
            idr_x1, idr_y1, idr_z1 = zip(*idr_points_1)
            ax.scatter(idr_x1, idr_y1, idr_z1, c=color1, s=20)

    if list2:
        x2, y2, z2 = zip(*list2)
        ax.plot(x2, y2, z2, color=color2, linewidth=1, label=name2)

        idr_indices_2 = [i for i, value in enumerate(IDR_2) if value > IDR_2_value]
        print(f"\n{name2}: atoms with IDR > {IDR_2_value}:")
        print(idr_indices_2)

        idr_points_2 = [list2[i] for i in idr_indices_2 if i < len(list2)]
        idr_points_2 = idr_points_2[:]

        if idr_points_2:
            idr_x2, idr_y2, idr_z2 = zip(*idr_points_2)
            ax.scatter(idr_x2, idr_y2, idr_z2, c=color2, s=20)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.legend(fontsize=8)
    ax.view_init(elev=30, azim=45, roll=15)
    plot.tight_layout()
    plot.show()
    #plot.savefig('figures/alignment.png')


def draw_aligned_records(records_file_path_1: str, color_1: str, name_1: str,
                         records_file_path_2: str, color_2: str, name_2: str):
    '''Draw a planar projection for a pair of aligned amyloid structures.'''

    points_1, seq_1, res_ids_1 = read_records(records_file_path_1)
    points_2, seq_2, res_ids_2 = read_records(records_file_path_2)

    seqres_1 = read_seqres(r'C:\Users\User\Documents\bioinf\smtb\seqres_all.txt', name_1)
    seqres_2 = read_seqres(r'C:\Users\User\Documents\bioinf\smtb\seqres_all.txt', name_2)

    aiupred_1 = read_aiupred(r'C:\Users\User\Documents\bioinf\smtb\AIUPred_all.txt', name_1)
    aiupred_2 = read_aiupred(r'C:\Users\User\Documents\bioinf\smtb\AIUPred_all.txt', name_2)

    IDR_1 = select_atomseq_idr(seqres_1, seq_1, aiupred_1)
    IDR_2 = select_atomseq_idr(seqres_2, seq_2, aiupred_2)

    print("\nIDR_1:")
    print(IDR_1)

    print("\nIDR_2:")
    print(IDR_2)

    print("\nLengths:")

    print("ATOM sequence 1:", len(seq_1))
    print("IDR_1:", len(IDR_1))

    print("ATOM sequence 2:", len(seq_2))
    print("IDR_2:", len(IDR_2))

    seq_1_aligned, seq_2_aligned = align_records(records_file_path_1,
                                                 records_file_path_2)

    paired_indices = get_aligned_indices(seq_1_aligned, seq_2_aligned)

    results = get_window_rmsd(seq_1_aligned, seq_2_aligned, points_1, points_2, WINDOW_SIZE)

    print("\nRMSD для рамок:")

    for frame, rmsd in results.items():
        print(f"{frame}: {rmsd}")

    mean_idr_values, rmsd_values = plot_rmsd_vs_mean_idr(results, IDR_1, IDR_2)

    #points_2 = best_rmsd_transform(points_1, points_2, paired_indices[-10:])

    #plane_basis = get_plane_basis(points_1, points_2, paired_indices[:])

    #_, xs_1, ys_1 = draw_projection(points_1, res_ids_1, seq_1, color_1, 10,
                                    #plane_basis)

    #_, xs_2, ys_2 = draw_projection(points_2, res_ids_2, seq_2, color_2, 0,
                                    #plane_basis)

    #connect_projected_points(xs_1, ys_1, xs_2, ys_2, paired_indices)

    #show_points_matplotlib(points_1, color_1, name_1, IDR_1, np.percentile(IDR_1, 50), points_2, color_2, name_2, IDR_2, np.percentile(IDR_2, 50))

    a = 0.3
    IDR_1_value = min(IDR_1[:]) + a * (max(IDR_1[:]) - min(IDR_1[:]))
    IDR_2_value = min(IDR_2[:]) + a * (max(IDR_2[:]) - min(IDR_2[:]))

    #show_points_matplotlib(points_1, color_1, name_1, IDR_1, IDR_1_value,
                           #points_2, color_2, name_2, IDR_2, IDR_2_value)


def read_structure_names(file_path: Path) -> list[str]:
    """Read structure names, one name per line."""
    print('?')
    with open(file_path, "r", encoding="utf8") as file:
        lines = file.readlines()

    structure_names = [line.strip() for line in lines[196:281]]
    
        
    print(structure_names)
    return structure_names

structure_names = read_structure_names(r"C:\Users\User\Documents\bioinf\smtb\families.txt")
structure_names.remove('9wap')

all_mean_idr = []
all_rmsd = []
WINDOW_SIZE = 10
OUTPUT_FILE = Path("figures/rmsd_vs_mean_idr_all_pairs.png")

for name_1, name_2 in combinations(structure_names, 2):

    print(name_1)
    print(name_2)

    extension_1 = ''
    extension_2 = ''

    folder = Path(r'C:\Users\User\Documents\bioinf\smtb\records')
    for file in folder.iterdir():
        if file.is_file() and file.name[:4] == name_1:
            extension_1 = file.name[4:8]
        if file.is_file() and file.name[:4] == name_2:
            extension_2 = file.name[4:8]

    records_file_path_1 = Path(r'C:\Users\User\Documents\bioinf\smtb\records\\' + name_1 + extension_1 + '_records.txt')
    records_file_path_2 = Path(r'C:\Users\User\Documents\bioinf\smtb\records\\' + name_2 + extension_2 + '_records.txt')


    with redirect_stdout(StringIO()):
        points_1, seq_1, res_ids_1 = read_records(records_file_path_1)
        points_2, seq_2, res_ids_2 = read_records(records_file_path_2)

        seqres_1 = read_seqres(r'C:\Users\User\Documents\bioinf\smtb\seqres_all.txt', name_1 + extension_1)
        seqres_2 = read_seqres(r'C:\Users\User\Documents\bioinf\smtb\seqres_all.txt', name_2 + extension_2)

        aiupred_1 = read_aiupred(r'C:\Users\User\Documents\bioinf\smtb\AIUPred_all.txt', name_1 + extension_1)
        aiupred_2 = read_aiupred(r'C:\Users\User\Documents\bioinf\smtb\AIUPred_all.txt', name_2 + extension_2)

        IDR_1 = select_atomseq_idr(seqres_1, seq_1, aiupred_1)
        IDR_2 = select_atomseq_idr(seqres_2, seq_2, aiupred_2)

        seq_1_aligned, seq_2_aligned = align_records(records_file_path_1,
                                                        records_file_path_2)
        
        paired_indices = get_aligned_indices(seq_1_aligned, seq_2_aligned)

        results = get_window_rmsd(seq_1_aligned, seq_2_aligned, points_1, points_2, WINDOW_SIZE)    

        mean_idr_values, rmsd_values = plot_rmsd_vs_mean_idr(results, IDR_1, IDR_2)

    all_mean_idr.extend(mean_idr_values)
    all_rmsd.extend(rmsd_values)

plot.figure()

plot.scatter(
    all_mean_idr,
    all_rmsd,
    s=40
)

plot.xlabel("Mean IDR")
plot.ylabel("RMSD")
plot.title(
    f"RMSD vs Mean IDR, window size = {WINDOW_SIZE}"
)

plot.tight_layout()

OUTPUT_FILE.parent.mkdir(
    parents=True,
    exist_ok=True
)

plot.savefig(
    OUTPUT_FILE,
    dpi=300
)

plot.close()

print("Completed")