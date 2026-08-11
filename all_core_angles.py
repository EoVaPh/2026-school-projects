from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from core_phi_psi_functions import (read_alignment, read_angles, collect_core_angles)


ALIGNMENTS_FOLDER = Path("MAFFT_families_atomseq_multiple_alignment")
ANGLES_FILE = Path("all_phi_psi.txt")
OUTPUT_FOLDER = Path("cos_sin_scatters")
OUTPUT_FOLDER.mkdir(exist_ok=True)


def save_scatter(data, ylabel, title, filename, core_length):
    """
    Create and save one scatter plot.

    Positive values  -> seagreen
    Negative values  -> steelblue
    Zero values      -> black
    """

    positions = np.arange(core_length)

    plt.figure(figsize=(16, 6))

    for position, values in enumerate(data):

        plt.axvline(
            position,
            color="lightgray",
            alpha=0.6,
            linewidth=0.8
        )
        positive = [value for value in values if value > 1e-10]
        negative = [value for value in values if value < -1e-10]
        zero = [value for value in values if abs(value) <= 1e-10]

        plt.scatter(
            [position] * len(negative),
            negative,
            color="steelblue",
            s=20
        )
        plt.scatter(
            [position] * len(positive),
            positive,
            color="seagreen",
            s=20
        )
        plt.scatter(
            [position] * len(zero),
            zero,
            color="black",
            s=20
        )

    plt.xlabel("Residue number in core")
    plt.ylabel(ylabel)
    plt.title(title)

    plt.xticks(
        positions,
        rotation=90
    )

    plt.ylim(-1.05, 1.05)

    plt.tight_layout()

    plt.savefig(
        filename,
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()


alignment_files = sorted(ALIGNMENTS_FOLDER.rglob("*_aligned.txt"))

for alignment_file in alignment_files:

    family_name = alignment_file.stem

    if family_name.endswith("_aligned"):
        family_name = family_name[:-8]


    print(f"Processing {family_name}...")

    family_output_folder = (OUTPUT_FOLDER / family_name)
    family_output_folder.mkdir(parents=True,exist_ok=True)

    alignment = read_alignment(alignment_file)

    if not alignment:
        print(
            f"  Alignment is empty: "
            f"{alignment_file}"
        )
        continue

    angles = read_angles(alignment, ANGLES_FILE)

    core, core_phis, core_psis = (collect_core_angles(alignment, angles, threshold=0.875))
    if not core:
        print(
            f"  Core was not found for "
            f"{family_name}"
        )
        continue
    print(core)


    core_cos_phis = [[np.cos(phi) for phi in phis if phi is not None] for phis in core_phis]
    core_cos_psis = [[np.cos(psi) for psi in psis if psi is not None] for psis in core_psis]
    core_sin_phis = [[np.sin(phi) for phi in phis if phi is not None] for phis in core_phis]
    core_sin_psis = [[np.sin(psi) for psi in psis if psi is not None] for psis in core_psis]

    save_scatter(
        core_cos_phis,
        "cos(φ)",
        f"{family_name}: cos(φ) angles in the core",
        family_output_folder
        / f"{family_name}_phi_cosine.png",
        len(core)
    )
    save_scatter(
        core_cos_psis,
        "cos(ψ)",
        f"{family_name}: cos(ψ) angles in the core",
        family_output_folder
        / f"{family_name}_psi_cosine.png",
        len(core)
    )
    save_scatter(
        core_sin_phis,
        "sin(φ)",
        f"{family_name}: sin(φ) angles in the core",
        family_output_folder
        / f"{family_name}_phi_sine.png",
        len(core)
    )
    save_scatter(
        core_sin_psis,
        "sin(ψ)",
        f"{family_name}: sin(ψ) angles in the core",
        family_output_folder
        / f"{family_name}_psi_sine.png",
        len(core)
    )

    print(
        f"  Saved 4 plots to "
        f"{family_output_folder}"
    )


print("Done.")