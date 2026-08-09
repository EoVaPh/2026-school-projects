from itertools import combinations
from pathlib import Path
import numpy as np
from contextlib import redirect_stdout
from io import StringIO

from draw_amyloid_function import best_rmsd_transform, read_records, get_aligned_indices, align_records, read_aiupred, read_seqres, select_atomseq_idr


def calculate_rmsd(points_1: list, points_2: list, pairs: list) -> list:
    transformed_points_2 = best_rmsd_transform(points_1, points_2, pairs)
    squared_distances = []

    for p1, p2 in zip(points_1, transformed_points_2):
        distance_squared = np.sum((np.array(p1) - np.array(p2)) ** 2)
        squared_distances.append(distance_squared)

    rmsd = np.sqrt(np.mean(squared_distances))

    return rmsd


def read_families(file_path: Path) -> dict:
    families = {}
    current_family = None

    with open(file_path, "r", encoding="utf8") as file:
        for line in file:
            line = line.strip()

            if not line:
                continue

            if line[0] != '>':
                families[current_family].append(line)
            else:
                current_family = line
                families[current_family] = []

    return families


families = read_families('families.txt')

file = open('all_rmsd.txt', 'w', encoding='utf8')

for family_name, structures in families.items():
    print()
    print(f"Family: {family_name}")
    print(f"Structures: {len(structures)}")

    #if len(structures) == 1:
    #    file.write('Only one structure')
    #else:
    if len(structures) > 1:
        pairs = list(combinations(structures, 2))

        file.write(f'{family_name}\n')

        count = 0

        for name_1, name_2 in pairs:
            count += 1
            print(f'{name_1}_{name_2}', count, '/', len(pairs))

            extension_1 = ''
            extension_2 = ''

            folder = Path('records')
            for record_file in folder.iterdir():
                if record_file.is_file() and record_file.name[:4] == name_1:
                    extension_1 = record_file.name[4:8]
                if record_file.is_file() and record_file.name[:4] == name_2:
                    extension_2 = record_file.name[4:8]

            records_file_path_1 = Path('records/' + name_1 + extension_1 + '_records.txt')
            records_file_path_2 = Path('records/' + name_2 + extension_2 + '_records.txt')

            with redirect_stdout(StringIO()):
                try:
                    points_1, seq_1, res_ids_1 = read_records(records_file_path_1)
                    points_2, seq_2, res_ids_2 = read_records(records_file_path_2)

                    seq_1_aligned, seq_2_aligned = align_records(records_file_path_1,
                                                                    records_file_path_2)

                    paired_indices = get_aligned_indices(seq_1_aligned, seq_2_aligned)

                    rmsd = calculate_rmsd(points_1, points_2, paired_indices)

                except Exception as error:
                    print(f"ERROR: {name_1} - {name_2}: {error}")

            file.write(f'{name_1}_{name_2}: {rmsd:.3f}\n')

file.close()

print('Completed')
