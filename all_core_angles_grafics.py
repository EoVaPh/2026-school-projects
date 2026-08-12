from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from core_phi_psi_functions import (read_alignment, read_angles, collect_core_angles, plot_scatter)


ALIGNMENTS_FOLDER = Path("MAFFT_families_atomseq_multiple_alignment")
ANGLES_FILE = Path("all_phi_psi.txt")
OUTPUT_FOLDER = Path("cos_sin_scatters")
OUTPUT_FOLDER.mkdir(exist_ok=True)


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

    core_cos_phis = [[np.cos(phi) for phi in phis if phi is not None] for phis in core_phis]
    core_cos_psis = [[np.cos(psi) for psi in psis if psi is not None] for psis in core_psis]
    core_sin_phis = [[np.sin(phi) for phi in phis if phi is not None] for phis in core_phis]
    core_sin_psis = [[np.sin(psi) for psi in psis if psi is not None] for psis in core_psis]

    positions = np.arange(len(core))

    plot_scatter(
        core_cos_phis,
        positions,
        "cos(φ)",
        f"{family_name}: cos(φ) angles in the core",
        family_output_folder
        / f"{family_name}_phi_cosine.png"
    )
    plot_scatter(
        core_cos_psis,
        positions,
        "cos(ψ)",
        f"{family_name}: cos(ψ) angles in the core",
        family_output_folder
        / f"{family_name}_psi_cosine.png"
    )
    plot_scatter(
        core_sin_phis,
        positions,
        "sin(φ)",
        f"{family_name}: sin(φ) angles in the core",
        family_output_folder
        / f"{family_name}_phi_sine.png"
    )
    plot_scatter(
        core_sin_psis,
        positions,
        "sin(ψ)",
        f"{family_name}: sin(ψ) angles in the core",
        family_output_folder
        / f"{family_name}_psi_sine.png"
    )

    print(
        f"  Saved 4 plots to "
        f"{family_output_folder}"
    )


print("Done.")