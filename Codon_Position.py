import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy.stats import chi2_contingency

#input files for the 4 mutation pools
mutation_files = {
    "iSNV_WH": "filtered_iSNV_WH_50x.csv",
    "iSNV_TOT": "filtered_iSNV_TOT_50x.csv",
    "SNP_WH": "filtered_SNP_WH_50x.csv",
    "SNP_TOT": "filtered_SNP_TOT_50x.csv"
}

mutation_classes = ["C>T", "G>A", "T>C", "A>G", "G>T", "A>T", "C>A", "G>C", "A>C", "T>G", "T>A"]

codon_position_counts = {pool: {1: 0, 2: 0, 3: 0} for pool in mutation_files.keys()}
mutation_position_counts = {pool: {mutation: {1: 0, 2: 0, 3: 0} for mutation in mutation_classes} for pool in mutation_files.keys()}

for pool, file_path in mutation_files.items():
    if not os.path.exists(file_path):
        print(f"Not found:{file_path}")
        continue

    df = pd.read_csv(file_path, sep=";")
    df["REF_CODON"] = df["REF_CODON"].astype(str).str.upper()
    df["ALT_CODON"] = df["ALT_CODON"].astype(str).str.upper()
    df = df[df["REF_CODON"].str.len() == 3]  # Keep only valid codons
    df = df[df["ALT_CODON"].str.len() == 3]

    # Identify the codon position that changed
    for _, row in df.iterrows():
        ref_codon, alt_codon = row["REF_CODON"], row["ALT_CODON"]
        for pos in range(3):
            if ref_codon[pos] != alt_codon[pos]:  # Detect mutation position
                codon_position_counts[pool][pos + 1] += 1  # Position 1-based
                mutation_type = f"{ref_codon[pos]}>{alt_codon[pos]}"  # Example: C>T
                if mutation_type in mutation_classes:
                    mutation_position_counts[pool][mutation_type][pos + 1] += 1
                break  # Only count the first detected change

# Chi-square test for overall mutation distribution across codon positions
chi2_results = {}
for pool in mutation_files.keys():
    observed = list(codon_position_counts[pool].values())
    expected = [sum(observed) / 3] * 3  # Assuming uniform distribution

    chi2_stat, p_value, _, _ = chi2_contingency([observed, expected])
    chi2_results[pool] = p_value  # Store only p-value

#first figure: Codon position mutation accumulation
fig, axes = plt.subplots(2, 2, figsize=(5, 5))
axes = axes.flatten()
colors = ["blue", "orange", "green", "red"]
titles = ["iSNV_WH", "iSNV_TOT", "SNP_WH", "SNP_TOT"]

for i, (pool, counts) in enumerate(codon_position_counts.items()):
    ax = axes[i]
    ax.bar(counts.keys(), counts.values(), color=colors[i], alpha=0.7)
    ax.set_title(f"({titles[i]})")
    ax.set_ylabel("Number of Mutations")
    ax.set_xticks([1, 2, 3])  # Codon positions


    # p-value result 
    p_value = chi2_results[pool]
    ax.text(2.7, max(counts.values()) * 0.9, f"p = {p_value:.4f}",
            ha="left", va="top", fontsize=8, bbox=dict(facecolor='white', alpha=0.5))

plt.tight_layout()
plt.show()

# second figure: Mutation classes per codon position
fig, axes = plt.subplots(4, 3, figsize=(15, 12)) 
axes = axes.flatten()

for i, mutation in enumerate(mutation_classes):
    ax = axes[i]
    width = 0.2  # Bar width for stacked bars

    for j, (pool, data) in enumerate(mutation_position_counts.items()):
        counts = [data[mutation][1], data[mutation][2], data[mutation][3]]  # Positions 1, 2, 3
        ax.bar([1 + j * width, 2 + j * width, 3 + j * width], counts, width=width, color=colors[j], alpha=0.7, label=titles[j])

    ax.set_title(mutation)
    ax.set_ylabel("Number of Mutations")
    ax.set_xticks([1.3, 2.3, 3.3])  
    ax.set_xticklabels(["1st", "2nd", "3rd"])
    ax.set_xlabel("Codon Position")

    if i == 0:  # Add legend to first plot only
        ax.legend(title="Mutation Pools")

plt.tight_layout()
plt.show()
