import numpy as np
from matplotlib import pyplot as plot
import matplotlib

from sequences_alignment import align_records
from alignment_index import get_aligned_indices

matplotlib.rcParams['figure.dpi'] = 300
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rc('font', family='STIXGeneral')
matplotlib.rc('font', weight='ultralight')


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


def connect_projected_points(xs_1: list, ys_1: list, xs_2: list, ys_2: list,
                             paired_indices: list):
    for pair in paired_indices:
        i, j = pair[0], pair[1]
        plot.plot([xs_1[i], xs_2[j]], [ys_1[i], ys_2[j]], linewidth=1,
                  color='black', zorder = -1, linestyle='--')


def show_points_matplotlib(list1, color1, name1, list2, color2, name2):
    fig = plot.figure()
    ax = fig.add_subplot(111, projection='3d')

    if list1:
        x1, y1, z1 = zip(*list1)
        ax.plot(x1, y1, z1, color=color1, linewidth=1, label=name1)
        #ax.scatter(x1, y1, z1, c=color1, s=20)
    if list2:
        x2, y2, z2 = zip(*list2)
        ax.plot(x2, y2, z2, color=color2, linewidth=1, label=name2)
        #ax.scatter(x2, y2, z2, c=color2, s=20)

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

    seq_1_aligned, seq_2_aligned = align_records(records_file_path_1,
                                                 records_file_path_2)

    paired_indices = get_aligned_indices(seq_1_aligned, seq_2_aligned)

    points_2 = best_rmsd_transform(points_1, points_2, paired_indices[:])

    #plane_basis = get_plane_basis(points_1, points_2, paired_indices[:])

    #_, xs_1, ys_1 = draw_projection(points_1, res_ids_1, seq_1, color_1, 10,
                                    #plane_basis)

    #_, xs_2, ys_2 = draw_projection(points_2, res_ids_2, seq_2, color_2, 0,
                                    #plane_basis)

    #connect_projected_points(xs_1, ys_1, xs_2, ys_2, paired_indices)

    show_points_matplotlib(points_1, color_1, name_1, points_2, color_2, name_2)


#draw_records('all_amyloid_structures/2beg.pdb_records.txt', '#f4acb7', 0)
#draw_records('all_amyloid_structures/2mxu.pdb_records.txt', '#669bbc', 1)

draw_aligned_records('all_amyloid_structures/9tpt.pdb_records.txt', '#f4acb7', '9tpt',
                     'all_amyloid_structures/7v49.pdb_records.txt', '#669bbc', '7v49')

plot.xticks([])
plot.yticks([])
plot.tight_layout()
plot.savefig('figures/proj.png')
