import numpy as np
from matplotlib import pyplot as plot
import matplotlib

matplotlib.rcParams['figure.dpi'] = 300
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rc("font", family="STIXGeneral")
matplotlib.rc("font", weight="ultralight")


def draw_projection(points: list, res_ids: list, color: str, zorder: int):
    pts = np.array(points, dtype=float)
    mean = np.mean(pts, axis=0)
    centered = pts - mean

    cov = np.cov(centered, rowvar=False)

    eigvals, eigvecs = np.linalg.eigh(cov)

    plane_basis = eigvecs[:, 1:]

    coords_2d = centered @ plane_basis

    projections = [tuple(coord) for coord in coords_2d]

    xs, ys = [coord[0] for coord in projections], [coord[1] for coord in projections]

    plot.scatter(xs, ys, color=color, s=128, zorder=zorder)

    for i in range(len(points) - 1):
        if res_ids[i+1] - res_ids[i] == 1:
            plot.plot([xs[i], xs[i+1]], [ys[i], ys[i+1]], color=color,
                      linewidth=4, zorder=zorder, alpha=0.5)


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


def draw_aligned_records(records_file_path_1: str, color_1: str,
                         records_file_path_2: str, color_2: str):
    '''Draw a planar projection for a pair of aligned amyloid structures.'''

    points_1, seq_1, res_ids_1 = read_records(records_file_path_1)
    points_2, seq_2, res_ids_2 = read_records(records_file_path_2)

    paired_indices = [
        (0, 6),
        (1, 7),
        (2, 8),
        (3, 9),
        (4, 10),
        (5, 11),
        (6, 12),
        (7, 13),
        (8, 14),
        (9, 15),
        (10, 16),
        (11, 17),
        (12, 18),
        (13, 19),
        (14, 20),
        (15, 21),
        (16, 22),
        (17, 23),
        (18, 24),
        (19, 25),
        (20, 26),
        (21, 27),
        (22, 28),
        (23, 29),
        (24, 30),
        (25, 31)
    ]

    points_2 = best_rmsd_transform(points_1, points_2, paired_indices)

    draw_projection(points_1, res_ids_1, color_1, zorder=0)
    draw_projection(points_2, res_ids_2, color_2, zorder=1)


#draw_records('all_amyloid_structures/2beg.pdb_records.txt', '#f4acb7', 0)
#draw_records('all_amyloid_structures/2mxu.pdb_records.txt', '#669bbc', 1)

draw_aligned_records('all_amyloid_structures/2beg.pdb_records.txt', '#f4acb7',
                     'all_amyloid_structures/2mxu.pdb_records.txt', '#669bbc')

plot.xticks([])
plot.yticks([])
plot.tight_layout()
plot.savefig('figures/proj.png')
