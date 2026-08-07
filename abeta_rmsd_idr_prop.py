from pathlib import Path
from itertools import combinations

import numpy as np
from matplotlib import pyplot as plot

from draw_amyloid import (
    read_records,
    read_seqres,
    read_aiupred,
    select_atomseq_idr,
    get_window_rmsd
)

from sequences_alignment import align_records


STRUCTURES_FILE = Path(
    r"C:\Users\User\Documents\bioinf\smtb\families.txt"
)

RECORDS_FOLDER = Path(
    r"C:\Users\User\Documents\bioinf\smtb\records"
)

SEQRES_FILE = Path(
    r"C:\Users\User\Documents\bioinf\smtb\seqres_all.txt"
)

AIUPRED_FILE = Path(
    r"C:\Users\User\Documents\bioinf\smtb\AIUPred_all.txt"
)

WINDOW_SIZE = 10


def read_structure_names(file_path):
    """Read structure names from a file."""

    with open(file_path, "r", encoding="utf8") as file:
            lines = file.readlines()
    
    structure_names = [line.strip() for line in lines[196:281]]

def find_records_file(structure_name):
    """Find records file for a structure."""

    files = list(
        RECORDS_FOLDER.glob(
            f"{structure_name}*_records.txt"
        )
    )

    if not files:
        raise FileNotFoundError(
            f"Records file not found: {structure_name}"
        )

    return files[0]


def read_structure(structure_name):
    """Read coordinates, sequence and IDR values."""

    records_file = find_records_file(structure_name)

    points, atom_seq, res_ids = read_records(
        str(records_file)
    )

    seqres = read_seqres(
        str(SEQRES_FILE),
        structure_name
    )

    aiupred = read_aiupred(
        str(AIUPRED_FILE),
        structure_name
    )

    idr = select_atomseq_idr(
        seqres,
        atom_seq,
        aiupred
    )

    return points, atom_seq, idr


def get_mean_idr(
    frame_start,
    frame_end,
    idr_1,
    idr_2
):
    """Calculate mean IDR for one frame."""

    values = []

    for i in range(frame_start, frame_end + 1):

        mean_pair_idr = (
            idr_1[i] + idr_2[i]
        ) / 2

        values.append(mean_pair_idr)

    return np.mean(values)


structure_names = read_structure_names(
    STRUCTURES_FILE
)

structures = {}

for name in structure_names:
    structures[name] = read_structure(name)


all_mean_idr = []
all_rmsd = []


for name_1, name_2 in combinations(
    structure_names,
    2
):

    points_1, seq_1, idr_1 = structures[name_1]
    points_2, seq_2, idr_2 = structures[name_2]

    records_file_1 = find_records_file(name_1)
    records_file_2 = find_records_file(name_2)

    seq_1_aligned, seq_2_aligned = align_records(
        str(records_file_1),
        str(records_file_2)
    )

    results = get_window_rmsd(
        seq_1_aligned,
        seq_2_aligned,
        points_1,
        points_2,
        WINDOW_SIZE
    )


    for frame_name, rmsd in results.items():

        if rmsd is None:
            continue

        start, end = map(
            int,
            frame_name.split(" - ")
        )

        mean_idr = get_mean_idr(
            start,
            end,
            idr_1,
            idr_2
        )

        all_mean_idr.append(mean_idr)
        all_rmsd.append(rmsd)


plot.figure(figsize=(8, 6))

plot.scatter(
    all_mean_idr,
    all_rmsd,
    s=40
)

plot.xlabel("Mean IDR")
plot.ylabel("RMSD")

plot.title(
    f"RMSD vs Mean IDR "
    f"(window size = {WINDOW_SIZE})"
)

plot.tight_layout()

plot.show()
