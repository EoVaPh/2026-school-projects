import numpy as np
from matplotlib import pyplot as plot
import matplotlib

matplotlib.rcParams["figure.dpi"] = 300
matplotlib.rcParams["mathtext.fontset"] = "stix"
matplotlib.rc("font", family="STIXGeneral")
matplotlib.rc("font", weight="ultralight")


def draw_records(records_file_path: str, color: str, zorder: int):
    records_file = open(records_file_path, 'r')
    records = records_file.readlines()
    records_file.close()

    points, res_ids = [], []

    for record in records:
        tokens = record.split()
        points.append((float(tokens[0]), float(tokens[1]), float(tokens[2])))
        res_ids.append(int(tokens[4]))

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
                      linewidth=4, zorder=zorder)


draw_records('all_amyloid_structures/2beg.pdb_records.txt', 'red', 0)
draw_records('all_amyloid_structures/2mxu.pdb_records.txt', 'blue', 1)

plot.xticks([])
plot.yticks([])
plot.tight_layout()
plot.savefig('figures/proj.png')
