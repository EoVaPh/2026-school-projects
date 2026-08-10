import os
import shutil
from pathlib import Path


def families_apd(pdb_ids: list, family_path: str):
    #os.mkdir(str("C:\Protein_Physics_SMTB\ProteinPhysics\APD\A\\" + family_path))
    #os.mkdir(str("C:\Protein_Physics_SMTB\ProteinPhysics\APD\contacts\\" + family_path))
    #os.mkdir(str("C:\Protein_Physics_SMTB\ProteinPhysics\APD\svg_logs\\" + family_path))

    structures = Path("C:\Protein_Physics_SMTB\ProteinPhysics\APD\structures")
    id1 = ""
    id2 = ""
    helix_1 = "python helix.py--chains A --twist -0.0 --rise 4.75 --axis_xy 0.0 0.0 --n_copies 2 --output "
    helix_2 = "python helix.py--chains A --twist -0.0 --rise 4.75 --axis_xy 0.0 0.0 --n_copies 2 --output "
    contacts_1 = "python contacts.py "
    contacts_2 = "python contacts.py "
    compare = "python compare.py--flip1 --flip2 --rot1 180 --rot2 70"
    id1_A_old = ""
    id1_A_new = ""
    id2_A_old = ""
    id2_A_new = ""
    id1_contacts_old = ""
    id1_contacts_old = ""
    id2_contacts_old = ""
    id2_contacts_old = ""
    log_old = ""
    svg_old = ""
    log_new = ""
    svg_old = ""
    with open("family_log.txt", "a", encoding="utf8") as t:
        print((family_path+"\n"), file=t)
    for i in pdb_ids:
        for j in pdb_ids:
            if i == j:
                break
            for item in structures.iterdir():
                item = str(item)
                if i in item:
                    id1 = item[-8:]
                if j in item:
                    id2 = item[-8:]
            print(id1, id2)
            id1_A_old = id1[0:-4] + "_A" + id1[-4:]
            id1_A_new = "A\\" + family_path + "\\" + id1_A_old
            id2_A_old = id2[0:-4] + "_A" + id2[-4:]
            id2_A_new = "A\\" + family_path + "\\" + id2_A_old
            id1_contacts_old = id1_A_old[-11:-4] + "_contacts.csv"
            id1_contacts_bild = id1_A_old[-11:-4] + "_contacts.bild"
            id1_contacts_new = "contacts\\" + family_path + "\\" + id1_contacts_old
            id2_contacts_old = id2_A_old[-11:-4] + "_contacts.csv"
            id2_contacts_bild = id2_A_old[-11:-4] + "_contacts.bild"
            id2_contacts_new = "contacts\\" + family_path + "\\" + id2_contacts_old
            helix_1 = helix_1[0:15] + " " + id1 + " " + helix_1[15:] + id1_A_old
            helix_2 = helix_2[0:15] + " " + id2 + " " + helix_2[15:] + id2_A_old
            contacts_1 += id1_A_old
            contacts_2 += id2_A_old
            compare = compare[0:17] + " " + id1_contacts_old + " " + id2_contacts_old + " " + compare[17:]
            log_old = id1_contacts_old[:-4] + "_vs_" + id2_contacts_old[:-4] + ".log"
            svg_old = id1_contacts_old[:-4] + "_vs_" + id2_contacts_old[:-4] + ".svg"
            log_new = "svg_logs\\" + family_path + "\\" + log_old
            svg_new = "svg_logs\\" + family_path + "\\" + svg_old

            id1_old = "C:\Protein_Physics_SMTB\ProteinPhysics\APD\structures\\" + id1
            id2_old = "C:\Protein_Physics_SMTB\ProteinPhysics\APD\structures\\" + id2
            shutil.copy(id1_old, "C:\Protein_Physics_SMTB\ProteinPhysics\APD")
            shutil.copy(id2_old, "C:\Protein_Physics_SMTB\ProteinPhysics\APD")

            folder_path = r"C:\Protein_Physics_SMTB\ProteinPhysics\APD"
            os.system(helix_1)
            os.system(helix_2)
            os.system(contacts_1)
            os.system(contacts_2)
            os.system(compare)

            family_log = []
            with (open(log_old, "r", encoding="utf8") as d):
                lines = [line.strip() for line in d.readlines()]
                residues = lines[4].split()[-1]
                right_left_flips = lines[25].split()[-1]
                apd = lines[31].split()[-1]
                apd_z = lines[32].split()[-1]
                pair_log = id1[:-4] + " " + id2[:-4] + " " + residues + " " + right_left_flips + " " + apd + " " + apd_z
            with open("family_log.txt", "a", encoding="utf8") as t:
                print((pair_log+"\n"), file = t)

            os.replace(id1_A_old, id1_A_new)
            os.replace(id2_A_old, id2_A_new)
            os.replace(id1_contacts_old, id1_contacts_new)
            os.replace(id2_contacts_old, id2_contacts_new)
            os.replace(log_old, log_new)
            os.replace(svg_old, svg_new)
            os.remove(id1_contacts_bild)
            os.remove(id2_contacts_bild)
            os.remove("C:\Protein_Physics_SMTB\ProteinPhysics\APD\\" + id1)
            os.remove("C:\Protein_Physics_SMTB\ProteinPhysics\APD\\" + id2)
            id1 = ""
            id2 = ""
            helix_1 = "python helix.py--chains A --twist -0.0 --rise 4.75 --axis_xy 0.0 0.0 --n_copies 2 --output "
            helix_2 = "python helix.py--chains A --twist -0.0 --rise 4.75 --axis_xy 0.0 0.0 --n_copies 2 --output "
            contacts_1 = "python contacts.py "
            contacts_2 = "python contacts.py "
            compare = "python compare.py--flip1 --flip2 --rot1 180 --rot2 70"
            id1_A_old = ""
            id1_A_new = ""
            id2_A_old = ""
            id2_A_new = ""
            id1_contacts_old = ""
            id1_contacts_old = ""
            id2_contacts_old = ""
            id2_contacts_old = ""
            log_old = ""
            svg_old = ""
            log_new = ""
            svg_old = ""
    return ("----------")

pdb_ids = []
with open("families.txt", "r", encoding="utf8") as f:
    lines = f.readlines()
    for _ in range(len(lines)):
        lines[_] = lines[_].replace("\n", "")
    family_path = lines[0].replace(">", "")
    for line in lines[1:]:
        if ">" in line:
            families_apd(pdb_ids, family_path)
            family_path = line[1:].replace(">", "")
            pdb_ids = []
        else:
            pdb_ids.append(line)