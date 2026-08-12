import numpy as np
from matplotlib import pyplot as plot
import matplotlib
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path

from draw_family_with_core_functions import *


def process_family(alignment_file, angles_file, records_folder, output_folder, threshold, core_size):
    """
    Process one protein family.

    The first structure in the alignment is used as the reference.
    Every other structure is aligned to the reference and saved
    as a separate PNG image.
    """

    # ---------------------------------------------------------
    # Read alignment
    # ---------------------------------------------------------

    alignment = read_alignment(alignment_file)

    # ---------------------------------------------------------
    # Get structure names from the alignment
    # ---------------------------------------------------------

    structure_names = get_family_structure_names(alignment_file)

    """ if len(structure_names) < 2:
        print(
            f"Skipping {alignment_file.name}: "
            f"less than 2 structures."
        )
        return """

    # ---------------------------------------------------------
    # Create output folder for this family
    # ---------------------------------------------------------

    family_name = alignment_file.stem

    family_output_folder = (
        output_folder / family_name
    )

    family_output_folder.mkdir(
        parents=True,
        exist_ok=True
    )

    # ---------------------------------------------------------
    # First structure = reference
    # ---------------------------------------------------------

    structure_name_1 = structure_names[0]

    input_file_1 = (
        records_folder /
        f"{structure_name_1}_records.txt"
    )

    if not input_file_1.exists():
        print(
            f"Records file not found: {input_file_1}"
        )
        return

    points_1, seq_1, res_ids_1 = read_records(
        input_file_1
    )

    # ---------------------------------------------------------
    # Process all other structures
    # ---------------------------------------------------------

    for structure_name_2 in structure_names:

        print()
        print("=" * 60)
        print(
            f"Family: {family_name}"
        )
        print(
            f"Reference: {structure_name_1}"
        )
        print(
            f"Structure: {structure_name_2}"
        )
        print("=" * 60)

        input_file_2 = (
            records_folder /
            f"{structure_name_2}_records.txt"
        )

        if not input_file_2.exists():
            print(
                f"Records file not found: {input_file_2}"
            )
            continue

        # -----------------------------------------------------
        # Read second structure
        # -----------------------------------------------------

        points_2, seq_2, res_ids_2 = read_records(
            input_file_2
        )

        # -----------------------------------------------------
        # Output file
        # -----------------------------------------------------

        save_path = (
            family_output_folder /
            f"{structure_name_2.replace('.pdb', '')}.png"
        )

        # -----------------------------------------------------
        # Align and draw
        # -----------------------------------------------------

        try:

            plot_pairwise_structures_with_core(
                structure_name_1=structure_name_1,
                structure_name_2=structure_name_2,
                alignment=alignment,
                angles_file=angles_file,
                points_1=points_1,
                points_2=points_2,
                res_ids_1=res_ids_1,
                res_ids_2=res_ids_2,
                threshold=threshold,
                structure_color_1="#669bbc",
                structure_color_2="#040404",
                core_size=core_size,
                save_path=save_path
            )

            print(
                f"Saved: {save_path}"
            )

        except Exception as e:

            print(
                f"Error processing "
                f"{structure_name_2}: {e}"
            )


alignment_folder = Path("MAFFT_families_atomseq_multiple_alignment")
records_folder = Path("records")
angles_file = Path("all_phi_psi.txt")
output_folder = Path("pairwise_results")


for alignment_file in sorted(alignment_folder.glob("*_aligned.txt")):

    process_family(
        alignment_file=alignment_file,
        angles_file=angles_file,
        records_folder=records_folder,
        output_folder=output_folder,
        threshold=0.875,
        core_size=80
    )