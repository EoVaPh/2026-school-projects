import ahocorasick

ag_automaton = ahocorasick.Automaton()
nag_automaton = ahocorasick.Automaton()

ag_motifs, non_ag_motifs = [], []
f = open('records/waltzdb_export.csv', 'r')
records = f.readlines()[1:]
f.close()

for record in records:
    tokens = record.split(',')
    motif, behav = tokens[0][1:-1], tokens[1][1:-1]

    if behav == 'non-amyloid':
        non_ag_motifs.append(motif)
    elif behav == 'amyloid':
        ag_motifs.append(motif)
    else:
        print('Each motif must be either "amyloid" or "non-amyloid".')
        exit(0)

for idx, key in enumerate(ag_motifs):
    ag_automaton.add_word(key, (idx, key))

ag_automaton.make_automaton()

for idx, key in enumerate(non_ag_motifs):
    nag_automaton.add_word(key, (idx, key))

nag_automaton.make_automaton()

# Read the sequences from FASTA file!
seqs = dict()
id_seq = ""
with open("records/seqres_all_new.txt", "r", encoding = "utf8") as f:
    for line in f:
        if ">" in line:
            id_seq = line[1:5]
        else:
            seqs[id_seq] = line.strip()

amyloid_count = 0
nonamyloid_count = 0

with open ("records/detected_motifs.txt", "w", encoding="utf8") as f:
    for id_seq, seq in seqs.items():
        for end_index, (insert_order, original_value) in ag_automaton.iter(seq):
            start_index = end_index - len(original_value) + 1
            print('amyloidogenic', start_index, end_index, original_value, id_seq, file=f)
            amyloid_count += 1
        for end_index, (insert_order, original_value) in nag_automaton.iter(seq):
            start_index = end_index - len(original_value) + 1
            print('non-amyloidogenic', start_index, end_index, original_value, id_seq, file=f)
            nonamyloid_count += 1
    print("amyloid_count:", amyloid_count)
    print("nonamyloid_count:", nonamyloid_count)