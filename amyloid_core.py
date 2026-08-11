import os
import biotite.structure as struc
from biotite.structure import lddt
import numpy as np

# ИЗМЕНЕНО: Функция теперь принимает относительные номера и сама находит поправку в файле
def parse_custom_coords_file(filepath, target_residues_relative):
    coords_dict = {}
    actual_residues = []
    if not os.path.exists(filepath):
        print(f"Warning: File {filepath} not found!")
        return [None] * len(target_residues_relative)

    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) >= 5:
                try:
                    x = float(parts[0])
                    y = float(parts[1])
                    z = float(parts[2])
                    res_num = int(parts[4])  # Column 5 is the actual residue number

                    coords_dict[res_num] = [x, y, z]
                except ValueError:
                    continue

    # Если файл пустой или не распарсился
    if not coords_dict:
        return [[None, None, None]] * len(target_residues_relative)

    # Находим минимальный номер остатка в файле (первый существующий остаток)
    first_res_in_file = min(coords_dict.keys())

    # Поправка: так как counter в выравнивании начинается с 1, реальный номер будет:
    # first_res_in_file + (relative_number - 1)
    offset = first_res_in_file - 1

    # Generate the final list of lists strictly in the order of core residues
    output_list = []
    for rel_num in target_residues_relative:
        actual_num = rel_num + offset  # Применяем поправку на стартовый номер
        actual_residues.append(actual_num)
        if actual_num in coords_dict:
            output_list.append(coords_dict[actual_num])
        else:
            output_list.append([None, None, None])
    #print(filepath, actual_residues)
    return output_list

# Function to calculating RMSD
def pair_rmsd(core_coord_1, core_coord_2):
    squared_distances = []
    for p1, p2 in zip(core_coord_1, core_coord_2):
        # Если одного из остатков нет в ядре, пропускаем его (опционально, зависит от вашей логики)
        if None in p1 or None in p2:
            continue
        distance_squared = np.sum((np.array(p1) - np.array(p2)) ** 2)
        squared_distances.append(distance_squared)

    if not squared_distances:
        return 0.0
    rmsd = np.sqrt(np.mean(squared_distances))
    return(rmsd)


# Function to calculating LDDT
def pair_lddt(core_coord_1, core_coord_2):

    length = len(core_coord_1)

    id1_obj = struc.AtomArray(length)
    id1_obj.coord = np.array(core_coord_1, dtype=float)
    id1_obj.coord.reshape(length, 3)

    id2_obj = struc.AtomArray(length)
    id2_obj.coord = np.array(core_coord_2, dtype=float)
    id2_obj.coord.reshape(length, 3)

    for obj in [id1_obj, id2_obj]:
        obj.res_id = np.arange(1, length + 1)
        obj.res_name = np.array(["GLY"] * length)
        obj.atom_name = np.array(["CA"] * length)
        obj.chain_id = np.array(["A"] * length)
        obj.hetero = np.array([False] * length)

    lddt_score = lddt(reference=id1_obj, subject=id2_obj, exclude_same_residue=False)
    return (lddt_score)


# Function to compare pair in core:
def compare_pair(dictionary_with_atoms):
    if len(dictionary_with_atoms) == 1:
        print("Family with only one member")
        with open("family_core_LDDT.txt", "a", encoding="utf8") as f:
            print("-family with only one member", file=f)
        with open("family_core_RMSD.txt", "a", encoding="utf8") as f:
            print("-family with only one member", file=f)
    for name_1 in dictionary_with_atoms:
        for name_2 in dictionary_with_atoms:
            if name_1 == name_2:
                break
            coords_1 = dictionary_with_atoms[name_1]
            coords_2 = dictionary_with_atoms[name_2]
            lddt_for_print = str(round(pair_lddt(coords_1, coords_2), 3))
            rmsd_for_print = str(round(pair_rmsd(coords_1, coords_2), 3))
            with open("family_core_LDDT.txt", "a", encoding="utf8") as f:
                print(name_1+" "+name_2+" "+lddt_for_print, file = f)
            with open("family_core_RMSD.txt", "a", encoding="utf8") as f:
                print(name_1+" "+name_2+" "+rmsd_for_print, file=f)
            with open("family_core_LDDT_RMSD.txt", "a", encoding="utf8") as f:
                print(name_1+" "+name_2+" "+" "+lddt_for_print+" "+rmsd_for_print, file = f)
            print(name_1, name_2, lddt_for_print, rmsd_for_print)

# Main loop: calculate core residue numbers and populate the dictionary
# Input alignment
folder_path = "C:\Protein_Physics_SMTB\ProteinPhysics\\amyloid_core\\families"  # путь к вашей папке
for filename in os.listdir(folder_path):
    print(filename)
    file_path = os.path.join(folder_path, filename)
    with open(file_path, "r", encoding="utf-8") as f:
        alignment_data = f.readlines()
    with open("family_core_LDDT.txt", "a", encoding="utf8") as f:
        print((">"+filename), file=f)
    with open("family_core_RMSD.txt", "a", encoding="utf8") as f:
        print((">"+filename), file=f)

    try:
        # Parse alignment lines
        sequences = {}
        for line in alignment_data:
            if line.strip():
                pdb_id, seq = line.split()
                sequences[pdb_id] = seq
        # Find core column indices in the alignment (where no sequence has a gap)
        sample_seqs = list(sequences.values())
        num_cols = len(sample_seqs[0])
        core_indices = [
            col for col in range(num_cols)
            if all(seq[col] != '-' for seq in sequences.values())
        ]
        if not core_indices:
            print("There is not a common core")
            with open("family_core_LDDT.txt", "a", encoding="utf8") as f:
                print(("-there is not a common core"), file=f)
            with open("family_core_RMSD.txt", "a", encoding="utf8") as f:
                print("-there is not a common core", file=f)
            continue

        start_col = min(core_indices)
        end_col = max(core_indices)

        all_structures_coords = {}
        for pdb_id, seq in sequences.items():
            target_residues_relative = []  # Относительные номера остатков (от 1 до длины структуры)
            amino_acid_counter = 0  # Track residue positions (1-based), ignoring gaps

            for col_idx, char in enumerate(seq):
                if char != '-':
                    amino_acid_counter += 1
                    # If the column is inside the alignment core, save the relative residue number
                    if start_col <= col_idx <= end_col:
                        target_residues_relative.append(amino_acid_counter)

            # Form file path: records/name_records.txt
            coords_filename = f"{pdb_id}_records.txt"
            full_path = os.path.join("C:\Protein_Physics_SMTB\ProteinPhysics\\amyloid_core\\records", coords_filename)
            #print(pdb_id, target_residues_relative)

            # Передаем относительные индексы, функция внутри пересчитает их с учетом реального старта
            all_structures_coords[pdb_id] = parse_custom_coords_file(full_path, target_residues_relative)

        # The all_structures_coords variable now contains the completed dictionary of lists
        print(all_structures_coords)
        compare_pair(all_structures_coords)
    except Exception as e:
        print(f"Failed: {e}")
        print(f"Type of error: {type(e).__name__}")
        with open("family_core_LDDT.txt", "a", encoding="utf8") as f:
            print((f"-failed: {e}"), file=f)
        with open("family_core_RMSD.txt", "a", encoding="utf8") as f:
            print((f"-failed: {e}"), file=f)
        continue