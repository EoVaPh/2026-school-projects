import os
import numpy as np
import matplotlib.pyplot as plt


def parse_family_file(file_path):
    """
    Parses the family file and extracts numerical pairwise comparison values.
    """
    families = {}
    current_family = None

    if not os.path.exists(file_path):
        print(f"Error: File {file_path} not found.")
        return {}

    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Detect the start of a new family
            if line.startswith('>'):
                current_family = line[1:].replace('_aligned.txt', '').strip()
                families[current_family] = []
            elif current_family:
                parts = line.split()
                if not parts:
                    continue
                try:
                    value = 1-float(parts[-1])
                    families[current_family].append(value)
                except ValueError:
                    continue

    # Keep ONLY families containing 10 or more pairwise comparisons
    clean_families = {k: v for k, v in families.items() if len(v) >= 10}
    return clean_families


def plot_protein_families(data, output_image='LDDT_plot.png'):
    """
    Plots floating boxes (Q1-Q3) with range whiskers (min/max), median, mean, and data points.
    """
    if not data:
        print("No data available to plot (all families might have been filtered out due to < 10 elements).")
        return

    # Synchronize family names and corresponding data lists
    fam_names = list(data.keys())
    data_list = [data[name] for name in fam_names]

    # Calculate statistical metrics
    means = [np.mean(d) for d in data_list]
    medians = [np.median(d) for d in data_list]
    mins = [np.min(d) for d in data_list]
    maxs = [np.max(d) for d in data_list]

    # Calculate 1st and 3rd quartiles for the floating box bounds
    q1_list = [np.percentile(d, 25) for d in data_list]
    q3_list = [np.percentile(d, 75) for d in data_list]
    box_heights = [q3 - q1 for q1, q3 in zip(q1_list, q3_list)]

    # Calculate asymmetric whiskers (min/max range relative to the mean value)
    lower_error = [m - mn for m, mn in zip(means, mins)]
    upper_error = [mx - m for m, mx in zip(means, maxs)]
    asymmetric_error = [lower_error, upper_error]

    x = np.arange(len(fam_names))

    # Adjust layout width dynamically based on the total number of families
    plt.figure(figsize=(max(10, len(fam_names) * 0.8), 6.5))

    # 1. Plot full range whiskers (min to max) centered on the mean values
    plt.errorbar(x, means, yerr=asymmetric_error, fmt='none', ecolor='#457b9d',
                 elinewidth=1.5, capsize=6, zorder=2, label='Full Range (min/max)')

    # 2. Plot floating boxes representing the Interquartile Range (from Q1 to Q3)
    # Using the bottom parameter allows the bar to float in mid-air
    plt.bar(x, box_heights, bottom=q1_list, color='#a8dadc', edgecolor='#457b9d',
            alpha=0.8, width=0.5, zorder=3, label='Interquartile Range (Q1 - Q3)')

    # 3. Plot median lines (horizontal dashes placed inside or near the boxes)
    plt.scatter(x, medians, color='#e63946', marker='_', s=400, linewidths=1.5,
                zorder=5, label='Median')

    # 4. Plot mean markers (green diamonds to distinguish them from the median)
    plt.scatter(x, means, color='#2a9d8f', marker='D', s=45,
                zorder=5, label='Mean Value')

    # 5. Plot individual pairwise lDDT scores with horizontal jitter
    for i, name in enumerate(fam_names):
        y_points = data[name]
        x_points = np.random.normal(i, 0.05, size=len(y_points))
        plt.scatter(x_points, y_points, color='#1d3557', alpha=0.4,
                    edgecolors='none', s=20, zorder=4,
                    label='Pairwise lDDT Scores' if i == 0 else "")

    # Dynamically generate X-axis labels with sample size (n=...) underneath the family name
    x_labels_with_n = [f"{name}\n(n={len(data[name])})" for name in fam_names]
    plt.xticks(x, x_labels_with_n, rotation=45, ha='right')

    # Chart styling and layout
    plt.ylabel('Similarity Score (lDDT)')
    plt.title('Distribution of Pairwise Comparisons by Protein Families (Filter: ≥ 10 elements)')

    # Strictly lock the lDDT range between 0.0 and 1.0
    plt.ylim(0.0, 1.0)

    plt.grid(axis='y', linestyle='--', alpha=0.5, zorder=1)
    plt.legend(loc='lower left', bbox_to_anchor=(0.0, 0.02))
    plt.tight_layout()

    # Save the resulting chart
    plt.savefig(output_image, dpi=150)
    plt.close()
    print(f"Chart successfully saved to: {output_image}")


if __name__ == '__main__':
    # Define your input data filename here
    input_file = 'family_core_LDDT.txt'

    if os.path.exists(input_file):
        data = parse_family_file(input_file)
        plot_protein_families(data)
    else:
        print(f"File {input_file} not found. Please specify the correct filename in the input_file variable.")