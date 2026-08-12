import numpy as np
from matplotlib import pyplot as plot
import matplotlib
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path

from core_phi_psi_functions import read_alignment, find_core, collect_core_angles, read_angles, calculate_cos_std
from draw_amyloid_function import best_rmsd_transform


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
    x, y, z = zip(*points)

    ax.plot(
        x,
        y,
        z,
        color=color,
        linewidth=1,
        label=name
    )

    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_zlabel("")

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])

    ax.legend(fontsize=8)

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

    structure_id = Path(structure_name)

    for name, aligned_sequence in alignment:

        alignment_id = Path(name)

        if alignment_id == structure_id:
            sequence = aligned_sequence
            break

    if sequence is None:
        raise ValueError(
            f"{structure_name} not found in alignment"
        )

    core_points = []
    core_residues = []
    core_alignment_positions = []

    residue_index = 0

    for alignment_position, aa in enumerate(sequence):

        if aa == "-":
            continue

        if alignment_position in core_positions:

            if residue_index < len(points):

                core_points.append(points[residue_index])
                core_residues.append(aa)
                core_alignment_positions.append(alignment_position)

        residue_index += 1

    return (core_points, core_residues, core_alignment_positions)


def get_core_pairs(alignment, structure_name_1, structure_name_2, core_positions):
    """
    Get corresponding residue index pairs for core alignment positions.
    """

    sequence_1 = None
    sequence_2 = None

    structure_id_1 = Path(structure_name_1)
    structure_id_2 = Path(structure_name_2)

    for name, sequence in alignment:

        alignment_id = Path(name)

        if alignment_id == structure_id_1:
            sequence_1 = sequence

        if alignment_id == structure_id_2:
            sequence_2 = sequence

    if sequence_1 is None:
        raise ValueError(f"{structure_name_1} not found in alignment")

    if sequence_2 is None:
        raise ValueError(f"{structure_name_2} not found in alignment")

    index_1 = 0
    index_2 = 0

    core_pairs = []

    for alignment_position, (aa1, aa2) in enumerate(zip(sequence_1, sequence_2)):

        if aa1 != "-" and aa2 != "-":

            if alignment_position in core_positions:
                core_pairs.append((index_1, index_2))

        if aa1 != "-":
            index_1 += 1

        if aa2 != "-":
            index_2 += 1

    return core_pairs


def split_points_by_residue_number(points, res_ids):
    """Split coordinates into continuous fragments.
    A new fragment starts when residue numbers are not consecutive.
    """

    fragments = []

    if not points:
        return fragments

    current_fragment = [points[0]]

    for i in range(1, len(points)):

        if res_ids[i] == res_ids[i - 1] + 1:
            current_fragment.append(points[i])

        else:
            fragments.append(current_fragment)
            current_fragment = [points[i]]

    fragments.append(current_fragment)

    return fragments


def plot_pairwise_structures_with_core(structure_name_1, structure_name_2, alignment, angles_file, points_1, points_2, res_ids_1, res_ids_2, threshold, structure_color_1, structure_color_2, core_size, save_path):
    """
    Draw two structures aligned using RMSD calculated only on the core.
    """

    core_positions = find_core(alignment, threshold)

    core_pairs = get_core_pairs(alignment, structure_name_1, structure_name_2, core_positions)

    if len(core_pairs) < 3:
        raise ValueError(
            "Not enough core pairs for RMSD transformation."
        )

    transformed_points_2 = best_rmsd_transform(points_1, points_2, core_pairs)
    fragments_2 = split_points_by_residue_number(transformed_points_2, res_ids_2)

    # Extract core points after transformation
    core_points_1 = [points_1[index_1] for index_1, index_2 in core_pairs]
    core_points_2 = [transformed_points_2[index_2] for index_1, index_2 in core_pairs]


    core_points_1 = np.asarray(core_points_1, dtype=float)
    core_points_2 = np.asarray(core_points_2, dtype=float)
    squared_distances = np.sum((core_points_1 - core_points_2) ** 2, axis=1)
    rmsd = np.sqrt(np.mean(squared_distances))

    # Convert coordinates to numpy arrays
    points_1 = np.asarray(points_1, dtype=float)

    transformed_points_2 = np.asarray(transformed_points_2, dtype=float)
    

    angles = read_angles(alignment, angles_file)
    core, core_phis, core_psis = (collect_core_angles(alignment, angles, threshold))
    std_cos_phi = calculate_cos_std(core_phis)

    core_points_1, core_residues_1, core_alignment_positions_1 = get_core_points(alignment, core_positions, structure_name_1, points_1)
    core_points_2, core_residues_2, core_alignment_positions_2 = get_core_points(alignment, core_positions, structure_name_2, transformed_points_2)

    # Plot
    fig = plot.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection="3d")

    # Structure 1
    """ ax.plot(
        points_1[:, 0],
        points_1[:, 1],
        points_1[:, 2],
        color=structure_color_1,
        linewidth=1,
        alpha=0.7
    ) """

    # Structure 2 after RMSD transformation
    for fragment in fragments_2:

        fragment = np.asarray(fragment, dtype=float)

        if len(fragment) < 2:
            continue

        ax.plot(
            fragment[:, 0],
            fragment[:, 1],
            fragment[:, 2],
            color=structure_color_2,
            linewidth=1,
            alpha=0.7
        )

    # Core points colored by STD of cos(phi)
    std_by_position = {position: std for position, std in zip(core_positions, std_cos_phi)}

    # Create one common color scale for both structures
    all_std_values = [
        std_by_position[position]
        for position in core_alignment_positions_1
        if position in std_by_position and not np.isnan(std_by_position[position])
    ]

    all_std_values += [
        std_by_position[position]
        for position in core_alignment_positions_2
        if position in std_by_position and not np.isnan(std_by_position[position])
    ]

    if all_std_values:

        norm = plot.Normalize(vmin=np.nanmin(all_std_values), vmax=np.nanmax(all_std_values))

        cmap_1 = LinearSegmentedColormap.from_list(
            "std_gradient_1",
            [
                "#2398ff",
                "#8ecae6",
                "#5743df"
            ]
        )

        cmap_2 = LinearSegmentedColormap.from_list(
            "std_gradient_2",
            [
                "#65c459",
                "#d0d95b",
                "#db6161"
            ]
        )

        # Structure 1 core
        core_points_1 = np.asarray(core_points_1, dtype=float)
        std_values_1 = np.asarray([std_by_position[position] for position in core_alignment_positions_1], dtype=float)

        colors_1 = cmap_1(norm(std_values_1))

        """ ax.scatter(
            core_points_1[:, 0],
            core_points_1[:, 1],
            core_points_1[:, 2],
            c=colors_1,
            s=core_size,
            alpha=0.7,
            depthshade=False
        ) """

        # Structure 2 core
        core_points_2 = np.asarray(core_points_2, dtype=float)
        std_values_2 = np.asarray([std_by_position[position] for position in core_alignment_positions_2], dtype=float)

        colors_2 = cmap_2(norm(std_values_2))

        ax.scatter(
            core_points_2[:, 0],
            core_points_2[:, 1],
            core_points_2[:, 2],
            c=colors_2,
            s=core_size,
            alpha=0.8,
            depthshade=False
        )

        for point, aa in zip(core_points_2, core_residues_2):

            x, y, z = point

            ax.text(
                x,
                y,
                z,
                aa,
                fontsize=6,
                fontweight="bold",
                ha="center",
                va="center"
            )


        # Residue labels
        """ for point, aa in zip(
            core_points_1,
            core_residues_1
        ):

            x, y, z = point

            ax.text(
                x,
                y,
                z,
                aa,
                fontsize=4,
                fontweight="normal"
            ) """


    # Colorbar
    """ sm = plot.cm.ScalarMappable(
        cmap=cmap_1,
        norm=norm
    )

    sm.set_array([])

    cbar = fig.colorbar(
        sm,
        ax=ax,
        pad=0.1
    )

    cbar.set_label("STD of cos(φ)") """

    sm = plot.cm.ScalarMappable(cmap=cmap_2, norm=norm)
    sm.set_array([])

    cbar_ax = fig.add_axes([0.78, 0.30, 0.025, 0.40])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label("STD of cos(φ)")

    # Legend
    """ structure_1_legend = Line2D(
        [0],
        [0],
        color=structure_color_1,
        linewidth=2,
        label=structure_name_1
    ) """

    structure_2_legend = Line2D(
        [0],
        [0],
        color=structure_color_2,
        linewidth=2,
        label=structure_name_2
    )

    ax.legend(
        handles=[
            #structure_1_legend,
            structure_2_legend,
        ],
        fontsize=8,
        loc="upper right",
        bbox_to_anchor=(0.85, 0.85)
    )


    ax.set_facecolor("white")

    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False

    ax.xaxis.pane.set_edgecolor("white")
    ax.yaxis.pane.set_edgecolor("white")
    ax.zaxis.pane.set_edgecolor("white")

    ax.grid(False)

    # Axes 
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])

    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_zlabel("")

    ax.xaxis.line.set_color("white")
    ax.yaxis.line.set_color("white")
    ax.zaxis.line.set_color("white")

    ax.view_init(
        elev=90,
        azim=0,
        roll=0
    )

    fig.subplots_adjust(left=0.02, right=0.88, bottom=0.02, top=0.98)

    fig.savefig(
        save_path,
        dpi=300,
        bbox_inches="tight",
        pad_inches=0
    )

    plot.close(fig)