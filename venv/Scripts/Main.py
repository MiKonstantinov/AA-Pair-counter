import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import webbrowser
from collections import Counter

# Your path to file
file_path = 'C:\\Users\\aligned_peptides.txt'  # Example 'C:\\aligned_peptides.txt'

# P2–P1 pair count
pair_counts = Counter()
p1_counts = Counter()
p2_counts = Counter()

with open(file_path, 'r') as f:
    for line in f:
        line = line.strip()
        if len(line) >= 6:
            p2 = line[4]
            p1 = line[5]
            if p1 != 'X' and p2 != 'X':
                pair_counts[(p2, p1)] += 1
                p1_counts[p1] += 1
                p2_counts[p2] += 1

# !You can change number of amino acids of P1 (p1_counts.most_common(*)) and P2 (p2_counts.most_common(*))
top_p1 = [aa for aa, _ in p1_counts.most_common(7)]
top_p2 = [aa for aa, _ in p2_counts.most_common(6)]

# DataFrame initialization
matrix = pd.DataFrame(index=top_p2, columns=top_p1, data=0)

# Filling table
for (p2, p1), count in pair_counts.items():
    if p2 in top_p2 and p1 in top_p1:
        matrix.at[p2, p1] = count

matrix = matrix.fillna(0).astype(int)

# Table print
print(f"{'P2':<3} {'P1':<3} {'Count':>6}")
for p2 in top_p2:
    for p1 in top_p1:
        count = matrix.at[p2, p1]
        if count > 0:
            print(f"{p2:<3} {p1:<3} {count:>6}")

# Heatmap building
plt.figure(figsize=(8, 6))
plt.rcParams['font.family'] = 'Arial'
sns.heatmap(matrix, annot=True, fmt="d", cmap="Reds", linewidths=0.5, cbar_kws={'label': 'Count'})
plt.xlabel("P1 (Cleavage site)")
plt.ylabel("P2 (Before cleavage site)")
plt.title("Top amino acid pairs (P2–P1)")
plt.tight_layout()
plt.savefig("top_amino_acid_pairs.png", dpi=300)
plt.savefig("top_amino_acid_pairs.svg", format="svg")
webbrowser.open("top_amino_acid_pairs.png")
