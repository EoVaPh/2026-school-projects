from pathlib import Path

from core_phi_psi_functions import read_alignment
from draw_family_with_core_functions import *

matplotlib.rcParams['figure.dpi'] = 300
matplotlib.rcParams['mathtext.fontset'] = 'stix'



alignment_file = Path("MAFFT_families_atomseq_multiple_alignment/" "family_004_aligned.txt")
angles_file = Path("all_phi_psi.txt")

structure_names = get_family_structure_names(alignment_file)
alignment = read_alignment(alignment_file)
structure_name_1 = structure_names[0]

input_file_1 = Path("records") / f"{structure_name_1}_records.txt"

points_1, seq_1, res_ids_1 = read_records(input_file_1)

output_folder = Path("family004_pairwise_structures")
output_folder.mkdir(exist_ok=True)

for structure_name_2 in structure_names:

    input_file_2 = Path("records") / f"{structure_name_2}_records.txt"

    if not input_file_2.exists():
        print(f"Records file not found: {input_file_2}")
        continue

    print()
    print("=" * 60)
    print(f"Processing: {structure_name_1} -> {structure_name_2}")
    print("=" * 60)

    points_2, seq_2, res_ids_2 = read_records(input_file_2)

    save_path = output_folder / (f"{structure_name_2.replace('.pdb', '')}_"
        f"aligned_to_{structure_name_1.replace('.pdb', '')}.png"
    )

    plot_pairwise_structures_with_core(
        structure_name_1=structure_name_1,
        structure_name_2=structure_name_2,
        alignment=alignment,
        angles_file=angles_file,
        points_1=points_1,
        points_2=points_2,
        res_ids_1=res_ids_1,
        res_ids_2=res_ids_2,
        threshold=0.875,
        structure_color_1="#669bbc",
        structure_color_2="#030303",
        core_size=80,
        save_path=save_path
    )

    print(f"Saved: {save_path}")