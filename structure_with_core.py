import numpy as np
from matplotlib import pyplot as plot
import matplotlib
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path

from core_phi_psi_functions import read_alignment, find_core, collect_core_angles, read_angles


matplotlib.rcParams['figure.dpi'] = 300
matplotlib.rcParams['mathtext.fontset'] = 'stix'



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


def get_family_structure_names(filename: Path) -> list:

    structure_names = []

    with open(filename, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            parts = line.split(maxsplit=1)

            if len(parts) != 2:
                continue

            structure_names.append(parts[0])

    return structure_names


def get_plane_basis(points: list):
    """
    Find the best-fit plane for one structure.
    """

    pts = np.array(points, dtype=float)
    mean = np.mean(pts, axis=0)
    centered = pts - mean
    cov = np.cov(centered, rowvar=False)
    eigvals, eigvecs = np.linalg.eigh(cov)
    normal = eigvecs[:, 0]
    plane_basis = eigvecs[:, 1:]

    return mean, normal, plane_basis


def draw_structure(points, name, color, sphere_size=20):
    """
    Draw one structure in 3D with spheres.
    """

    points = np.asarray(points, dtype=float)

    fig = plot.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Structure
    x, y, z = zip(*points)

    ax.plot(
        x,
        y,
        z,
        color=color,
        linewidth=1,
        label=name
    )

    # Spheres
    #ax.scatter(x, y, z, color=color, s=sphere_size, depthshade=True)

    # Remove axes labels and ticks
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_zlabel("")

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])

    ax.legend(fontsize=8)

    # Same view as original code
    ax.view_init(
        elev=30,
        azim=45,
        roll=15
    )

    plot.tight_layout()
    plot.show()


def get_core_points(alignment, core_positions, structure_name, points):
    """
    Get coordinates of core residues for one structure.
    """

    sequence = None

    for name, aligned_sequence in alignment:

        if name == structure_name:
            sequence = aligned_sequence
            break

    if sequence is None:
        raise ValueError(
            f"{structure_name} not found in alignment"
        )

    core_points = []
    core_residues = []
    core_alignment_positions = []

    # Index in the original structure.
    # It increases only when we encounter a real residue.
    residue_index = 0

    for alignment_position, aa in enumerate(sequence):

        # Gap in THIS structure.
        # It has no coordinate and no amino-acid label.
        if aa == "-":
            continue

        # This position belongs to the family-wide core
        if alignment_position in core_positions:

            if residue_index < len(points):

                core_points.append(
                    points[residue_index]
                )

                # IMPORTANT:
                # Take the amino acid from THIS structure's
                # alignment sequence.
                core_residues.append(aa)

                # Save the alignment position so that
                # phi/psi statistics can be matched correctly.
                core_alignment_positions.append(
                    alignment_position
                )

        residue_index += 1

    return (
        core_points,
        core_residues,
        core_alignment_positions
    )


def calculate_cos_std(phi_angles):
    
    std_values = []

    for angles in phi_angles:
        angles = np.asarray(angles, dtype=float)
        angles = angles[~np.isnan(angles)]

        if len(angles) == 0:
            std_values.append(np.nan)
            continue

        cos_values = np.cos(angles)
        std = np.std(cos_values)
        std_values.append(std)

    return std_values


def plot_structure_with_core(structure_name, alignment, points, threshold, structure_color, core_size):
    """
    Draw one 3D protein structure.

    The whole structure is shown as a line.
    Core residues are highlighted with scatter points.
    """

    core_positions = find_core(alignment, threshold)
    core_points, core_residues, core_alignment_positions = get_core_points(alignment, core_positions, structure_name, points)
    angles = read_angles(alignment, angles_file)
    core, core_phis, core_psis = (collect_core_angles(alignment, angles, threshold))
    print(core)
    std_cos_phi = calculate_cos_std(core_phis)

    points = np.asarray(points, dtype=float)

    fig = plot.figure(figsize=(8, 8))

    ax = fig.add_subplot(111, projection="3d")

    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2]

    ax.plot(
        x,
        y,
        z,
        color=structure_color,
        linewidth=1,
        alpha=0.7
    )


    if core_points:

        core_points = np.asarray(core_points, dtype=float)
        std_by_position = {position: std for position, std in zip(core_positions, std_cos_phi)}
        std_values = np.asarray([std_by_position[position] for position in core_alignment_positions], dtype=float)
        norm = plot.Normalize(vmin=np.nanmin(std_values), vmax=np.nanmax(std_values))
        cmap = cmap = LinearSegmentedColormap.from_list("std_gradient", ["#65c459", "#d0d95b", "#db6161"])
        colors = cmap(norm(std_values))

        ax.scatter(
            core_points[:, 0],
            core_points[:, 1],
            core_points[:, 2],
            c=colors,
            s=core_size,
            depthshade=False
        )

        for point, aa in zip(core_points, core_residues):

            x, y, z = point

            ax.text(
                x,
                y,
                z,
                aa,
                fontsize=4,
                fontweight="normal"
            )
    
        sm = plot.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, pad=0.1)
        cbar.set_label("STD of cos(φ)")


    structure_legend = Line2D([0],
        [0],
        color=structure_color,
        linewidth=2,
        label=structure_name
    )

    core_legend = Line2D(
        [0],
        [0],
        marker="o",
        color="none",
        markerfacecolor=structure_color,
        markersize=6,
        label="Core"
    )

    ax.legend(
        handles=[
            structure_legend,
            core_legend
        ],
        fontsize=8
    )

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])

    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_zlabel("")

    ax.view_init(
        elev=30,
        azim=45,
        roll=15
    )

    plot.tight_layout()

    plot.show()



#structure_names = get_family_structure_names(alignment_file)

alignment_file = Path("MAFFT_families_atomseq_multiple_alignment\\family_004_aligned.txt")
angles_file = Path("all_phi_psi.txt")
input_file = "records/9tke.pdb_records.txt"

points, seq, res_ids = read_records(input_file)
alignment = read_alignment(alignment_file)


plot_structure_with_core(
    structure_name="9tke.pdb",
    alignment=alignment,
    points=points,
    threshold=0.875,
    structure_color= '#669bbc',
    core_size=20
)