import os
import biotite.structure as struc
from biotite.structure import lddt
import numpy as np

def get_paired_indices_core(seq_1: str, seq_2: str, core_begin: int, core_end: int) -> list:
    N = len(seq_1)

    assert len(seq_2) == N

    i, i_1, i_2 = 0, -1, -1
    paired_indices_core = []

    for i in range(N):
        if seq_1[i] != '-':
            i_1 += 1

        if seq_2[i] != '-':
            i_2 += 1

        if seq_1[i] != '-' and seq_2[i] != '-' and core_begin <= i <= core_end:
            paired_indices_core.append((i_1, i_2))

    return paired_indices_core

def best_rmsd_transform(list1: list, list2: list) -> list:
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

    Returns
    -------
    transformed : list of tuple of (float, float, float)
        All points from list2 after applying the optimal rotation and
        translation.
    """
    if list1 == [] or list2 == []:
        return None
    # Convert to numpy arrays (Nx3)
    ref = np.array(list1, dtype=float)
    mob = np.array(list2, dtype=float)
    pairs = []
    for i in range(len(list1)):
        pair = [i, i]
        pairs.append(pair)
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

def pair_rmsd(core_coord_1, core_coord_2):
    if core_coord_1 == [] or core_coord_2 == []:
        return "None"

    squared_distances = []
    core_coord_2 = best_rmsd_transform(core_coord_1, core_coord_2)

    for p1, p2 in zip(core_coord_1, core_coord_2):
        distance_squared = np.sum((np.array(p1) - np.array(p2)) ** 2)
        squared_distances.append(distance_squared)

    if not squared_distances:
        return "None"
    rmsd = np.sqrt(np.mean(squared_distances))
    return rmsd

def pair_lddt(core_coord_1, core_coord_2):

    length = len(core_coord_1)
    if length == 0:
        return 0.0

    id1_obj = struc.AtomArray(length)
    id1_obj.coord = np.array(core_coord_1, dtype=float).reshape((length, 3))

    id2_obj = struc.AtomArray(length)
    id2_obj.coord = np.array(core_coord_2, dtype=float).reshape((length, 3))

    for obj in [id1_obj, id2_obj]:
        obj.res_id = np.arange(1, length + 1)
        obj.res_name = np.array(["GLY"] * length)
        obj.atom_name = np.array(["CA"] * length)
        obj.chain_id = np.array(["A"] * length)
        obj.hetero = np.array([False] * length)

    try:
        lddt_score = lddt(reference=id1_obj, subject=id2_obj, exclude_same_residue=False)
        return float(np.mean(lddt_score))
    except Exception:
        return 0.0

def find_core(alignment, threshold=0.875):
    if not alignment:
        return []

    alignment_length = len(alignment[0][1])
    number_of_sequences = len(alignment)

    best_start = None
    best_end = None
    best_length = 0
    current_start = None
    current_length = 0

    for position in range(alignment_length):
        number_of_residues = sum(
            sequence[1][position] != "-" for sequence in alignment
        )
        fraction = number_of_residues / number_of_sequences

        if fraction >= threshold:
            if current_start is None:
                current_start = position
                current_length = 1
            else:
                current_length += 1

            if current_length > best_length:
                best_length = current_length
                best_start = current_start
                best_end = position
        else:
            current_start = None
            current_length = 0

    if best_start is None:
        return []
    return list(range(best_start, best_end + 1))


# Main loop
folder_path = "C:\\Protein_Physics_SMTB\\ProteinPhysics\\amyloid_core\\families"

for out_f in ["family_core_LDDT.txt", "family_core_RMSD.txt", "family_core_LDDT_RMSD.txt"]:
    with open(out_f, "w", encoding="utf8") as f:
        pass

for filename in os.listdir(folder_path):
    print(">"+filename[:-4])
    for out_f in ["family_core_LDDT.txt", "family_core_RMSD.txt", "family_core_LDDT_RMSD.txt"]:
        with open(out_f, "a", encoding="utf8") as f:
            print((">" + filename[:-4]), file=f)
    file_path = os.path.join(folder_path, filename)
    with open(file_path, "r", encoding="utf-8") as f:
        alignment_data = f.readlines()
    sequences = {}
    for line in alignment_data:
        if line.strip():
            pdb_id, seq = line.split()
            sequences[pdb_id] = seq
    if len(alignment_data) == 1:
        for out_f in ["family_core_LDDT.txt", "family_core_RMSD.txt", "family_core_LDDT_RMSD.txt"]:
            with open(out_f, "a", encoding="utf8") as f:
                print("-there is only one member of family", file=f)
    alignment_list = list(sequences.items())
    core_indices = find_core(alignment_list, threshold=0.875)

    if not core_indices:
        print("-there is not a common core")
        for out_f in ["family_core_LDDT.txt", "family_core_RMSD.txt", "family_core_LDDT_RMSD.txt"]:
            with open(out_f, "a", encoding="utf8") as f:
                print("-there is not a common core", file=f)
        continue

    start_col = min(core_indices)
    end_col = max(core_indices)
    print(f"Core: {start_col, end_col}")

    for id1 in sequences.keys():
        for id2 in sequences.keys():
            if id1 == id2:
                break
            pair_indeces = get_paired_indices_core(sequences[id1], sequences[id2], start_col, end_col)
            id1_indeces = []
            id2_indeces = []
            for i in range(len(pair_indeces)):
                id1_indeces.append(pair_indeces[i][0])
                id2_indeces.append(pair_indeces[i][1])
            with open(f"C:\Protein_Physics_SMTB\ProteinPhysics\\amyloid_core\\records\\{id1}_records.txt", "r", encoding="utf8") as f:
                atomseq = f.readlines()
                id1_atoms_coords = []
                for res_num in id1_indeces:
                    res = atomseq[res_num].split()
                    coords = []
                    x_coord = float(res[0])
                    coords.append(x_coord)
                    y_coord = float(res[1])
                    coords.append(y_coord)
                    z_coord = float(res[2])
                    coords.append(z_coord)
                    id1_atoms_coords.append(coords)
            with open(f"C:\Protein_Physics_SMTB\ProteinPhysics\\amyloid_core\\records\\{id2}_records.txt", "r", encoding="utf8") as f:
                atomseq = f.readlines()
                id2_atoms_coords = []
                for res_num in id2_indeces:
                    res = atomseq[res_num].split()
                    coords = []
                    x_coord = float(res[0])
                    coords.append(x_coord)
                    y_coord = float(res[1])
                    coords.append(y_coord)
                    z_coord = float(res[2])
                    coords.append(z_coord)
                    id2_atoms_coords.append(coords)
            with open("family_core_LDDT.txt", "a", encoding="utf8") as f:
                LDDT_score = str(round(pair_lddt(id1_atoms_coords, id2_atoms_coords), 3))
                lddt_for_print = id1+" "+id2+" "+LDDT_score
                print(lddt_for_print, file=f)
            with open("family_core_RMSD.txt", "a", encoding="utf8") as f:
                RMSD_score = pair_rmsd(id1_atoms_coords, id2_atoms_coords)
                if isinstance(RMSD_score, float):
                    RMSD_score = str(round(RMSD_score, 1))
                rmsd_for_print = id1+" "+id2+" "+RMSD_score
                print(rmsd_for_print, file=f)
            with open("family_core_LDDT_RMSD.txt", "a", encoding="utf8") as f:
                lddt_rmsd_for_print = id1+" "+id2+" "+LDDT_score+" "+RMSD_score
                print(lddt_rmsd_for_print, file=f)