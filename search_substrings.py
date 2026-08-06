import ahocorasick

ag_automaton = ahocorasick.Automaton()
nag_automaton = ahocorasick.Automaton()

ag_motifs, non_ag_motifs = [], []
f = open('waltzdb_export.csv', 'r')
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
        raise ValueException('Each motif must be either "amyloid" or "non-amyloid".')

for idx, key in enumerate(ag_motifs):
    ag_automaton.add_word(key, (idx, key))

ag_automaton.make_automaton()

for idx, key in enumerate(non_ag_motifs):
    nag_automaton.add_word(key, (idx, key))

nag_automaton.make_automaton()

# Read the sequences from FASTA file!
for seq in seqs:
    for end_index, (insert_order, original_value) in ag_automaton.iter(seq):
        print('amyloidogenic', end_index, insert_order, original_value)
    for end_index, (insert_order, original_value) in nag_automaton.iter(seq):
        print('non-amyloidogenic', end_index, insert_order, original_value)
