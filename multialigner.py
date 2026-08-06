from pathlib import Path
import seqres_alignment_algorithm
import networkx as nx
import time
from tqdm import tqdm

folder_path = Path("C:\\PP\\Protein Physics\\records")
for path1 in tqdm(folder_path.glob("*.fa"), desc = "progress"):
    for path2 in folder_path.glob("*.fa"):
        if path1 == path2:
            break
        else:
            seqres_alignment_algorithm.aligners(path1, path2)

data = seqres_alignment_algorithm.bypair

line = ""
with open("records/data.txt", "w", encoding="utf-8") as f:
    for i in data:
        line += str(i[0])+ " " + str(i[1]) + " " + str(i[2]) + "\n"
        f.write(line)
        line = ""