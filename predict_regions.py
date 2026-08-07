from aiupred import AIUPred
import numpy as np

# Define the target FASTA file.
target = 'C:\\Users\\User\\Documents\\bioinf\\smtb\\seqres_all_new.txt'

# Initialize the predictor.
predictor = AIUPred()

# Read all lines from the FASTA file.
fasta_file = open(target, 'r')
lines = fasta_file.readlines()
fasta_file.close()

# Number of sequences in the FASTA file.
num_seqs = len(lines) / 2

# Start writing to file containing the same annotations but pre-calculated IDR
# propensities instead of sequences.
output_fasta_file = open(r'C:\Users\User\Documents\bioinf\smtb\AIUPred_new.txt', 'w')

count = 0

for line in lines:
    if '>' in line:
        output_fasta_file.write(line)
    else:
        sequence = line.strip()

        disorder_propensities = predictor.predict_disorder(sequence)

        arr_str = np.array2string(disorder_propensities, threshold=np.inf)
        # Remove line separations and array brackets, create a single line
        # ending with new line symbol.
        output_fasta_file.write(
            arr_str.replace('\n', '').replace('[', '').replace(']', '') + '\n'
        )

        count += 1
        print(count, '/', num_seqs)

output_fasta_file.close()
