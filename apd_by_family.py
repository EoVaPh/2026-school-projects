import os
from pathlib import Path

with open ("test.txt", "r", encoding = "utf8") as f:
    line = f.readline().replace("\n", "").replace(">","")
    #os.mkdir(str("C:\Protein_Physics_SMTB\ProteinPhysics\APD\A\\")+line)
    #os.mkdir(str("C:\Protein_Physics_SMTB\ProteinPhysics\APD\contacts\\")+line)
    #os.mkdir(str("C:\Protein_Physics_SMTB\ProteinPhysics\APD\svg_logs\\")+line)
    pdb_ids = f.readlines()

for i in range(len(pdb_ids)):
    pdb_ids[i] = pdb_ids[i].replace("\n", "")

pdb_ids = ["9b4l","9b4m"]
structures = Path("C:\Protein_Physics_SMTB\ProteinPhysics\APD\structures")
id1 = ""
id2 = ""
helix_1 = "python helix.py--chains A --twist -1 --rise 4.75 --axis_xy 0.0 0.0 --n_copies 2 --output "
helix_2 = "python helix.py--chains A --twist -1 --rise 4.75 --axis_xy 0.0 0.0 --n_copies 2 --output "
contacts_1 = "python contacts.py "
contacts_2 = "python contacts.py "
compare = "python compare.py--flip1 --flip2 --rot1 180 --rot2 70"

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
        id1_A_old = id1[0:-4]+"_A"+id1[-4:]
        id1_A_new = "A\FILAMENT OF TAU IN COMPLEX WITH D-TLKIVWI, A D-PEPTIDE THAT\\" + id1_A_old
        id2_A_old = id2[0:-4] + "_A" + id2[-4:]
        id2_A_new = "A\FILAMENT OF TAU IN COMPLEX WITH D-TLKIVWI, A D-PEPTIDE THAT\\" + id2_A_old
        id1_contacts_old = id1_A_old[-11:-4] + "_contacts.csv"
        id1_contacts_new = "contacts\FILAMENT OF TAU IN COMPLEX WITH D-TLKIVWI, A D-PEPTIDE THAT\\" + id1_contacts_old
        id2_contacts_old = id2_A_old[-11:-4] + "_contacts.csv"
        id2_contacts_new= "contacts\FILAMENT OF TAU IN COMPLEX WITH D-TLKIVWI, A D-PEPTIDE THAT\\" + id2_contacts_old
        helix_1 = helix_1[0:15] + " " + id1 + " " + helix_1[15:] + id1_A_old
        helix_2 = helix_2[0:15] + " " + id2 + " " +  helix_2[15:] + id2_A_old
        contacts_1 += id1_A_old
        contacts_2 += id2_A_old
        compare = compare[0:17] + " " + id1_contacts_old + " " + id2_contacts_old + " " + compare[17:]
        print(helix_1)
        folder_path = r"C:\\Protein_Physics_SMTB\\ProteinPhysics\\APD"
        os.system(helix_1)
        os.system(helix_2)
        os.system(contacts_1)
        os.system(contacts_2)
        os.system(compare)
        id1_A_old = ""
        id1_A_new = ""
        id2_A_old = ""
        id2_A_new = ""
        id1_contacts_old = ""
        id1_contacts_old = ""
        id2_contacts_old = ""
        id2_contacts_old = ""
        helix_1 = "python helix.py--chains A --twist -1 --rise 4.75 --axis_xy 0.0 --n_copies 2 --output "
        helix_2 = "python helix.py--chains A --twist -1 --rise 4.75 --axis_xy 0.0 --n_copies 2 --output "
        contacts_1 = "python contacts.py "
        contacts_2 = "python contacts.py "
        compare = "python compare.py--flip1 --flip2 --rot1 180 --rot2 70"