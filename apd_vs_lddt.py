from pathlib import Path

from matplotlib import pyplot as plot
import matplotlib

matplotlib.rcParams['figure.dpi'] = 300
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rc('font', family='STIXGeneral')
matplotlib.rc('font', weight='ultralight')

from scipy import stats
import numpy as np


def pearson_r_ci(r, n, confidence=0.95):
    """
    Calculates the confidence interval for a Pearson r correlation coefficient.

    Parameters:
    r (float): The calculated correlation coefficient (r-value)
    n (int): The sample size
    confidence (float): The confidence level (default 0.95)
    """
    # 1. Fisher's z-transformation of the r-value
    z_prime = np.arctanh(r)

    # 2. Calculate the standard error of the transformed z
    stderr = 1 / np.sqrt(n - 3)

    # 3. Find the critical z-score for the 95% confidence level
    alpha = 1 - confidence
    z_critical = stats.norm.ppf(1 - alpha / 2)

    # 4. Compute the lower and upper bounds in z-space
    z_lower = z_prime - (z_critical * stderr)
    z_upper = z_prime + (z_critical * stderr)

    # 5. Inverse Fisher's z-transformation to get bounds back to r-space
    r_lower = np.tanh(z_lower)
    r_upper = np.tanh(z_upper)

    return r_lower, r_upper


family_core_LDDT_file = open('family_core_LDDT.txt', 'r')
family_core_LDDT_lines = family_core_LDDT_file.readlines()
family_core_LDDT_file.close()

LDDTs = dict()

for line in family_core_LDDT_lines:
    if '>' not in line and '-' not in line:
        tokens = line.split()
        min_id, max_id = min(tokens[0], tokens[1]), max(tokens[0], tokens[1])
        LDDTs[(min_id[:-4], max_id[:-4])] = float(tokens[2])

apd_path = Path('APD_output')

apd_files = set([f.name for f in apd_path.iterdir() if f.is_file()])

data = []

for apd_file in apd_files:
    if '_vs_' in apd_file and '.svg' in apd_file:
        apd_file_tokens = apd_file.split('_')
        id_1, id_2 = apd_file_tokens[1], apd_file_tokens[5]

        apd_file_lines = open('APD_output/' + apd_file, 'r').readlines()
        for line in apd_file_lines:
            if 'residues in common:' in line:
                common_residues = int(line.strip().replace('</tspan></text>', '').split('>')[-1])

            if 'right-left flips:' in line:
                flips = int(line.strip().replace('</tspan></text>', '').split('>')[-1])

                min_id, max_id = min(id_1, id_2), max(id_1, id_2)

                if (min_id, max_id) in LDDTs:
                    data.append((flips / common_residues, LDDTs[(min_id, max_id)]))

                break


flip_data = [entry[0] for entry in data]
LDDT_data = [1 - entry[1] for entry in data]

plot.scatter(flip_data, LDDT_data, color='#e36414', alpha=0.5)

slope, intercept, r_value, p_value, std_err = stats.linregress(flip_data,
                                                               LDDT_data)

arange = np.arange(min(flip_data), max(flip_data), 0.001)
plot.plot(arange, arange * slope + intercept, color='#a44200')

CI = pearson_r_ci(r_value, len(flip_data))

plot.ylim(ymax=0.7)

plot.xlabel('(Number of right-left flips) / (Number of common residues)', fontsize=16)
plot.ylabel('(1 – LDDT) of core substructures', fontsize=16)
plot.title(r'r$^2$ = ' + str(round(r_value**2, 3)) + ', CI = [' + str(round(CI[0], 3)) + ', ' + str(round(CI[1], 3)) + ']')
plot.tight_layout()
plot.savefig('APD_vs_LDDT.png')
