from Bio import AlignIO
from Bio.Seq import Seq
from collections import Counter
import matplotlib.pyplot as plt
import pandas as pd
import os

# Load nucleotide alignment
alignment = AlignIO.read(
    "data/sequences/HIV1_protease_aligned.fasta",
    "fasta"
)

protein_sequences = []

# Translate aligned nucleotide sequences
for record in alignment:
    nucleotide = str(record.seq).upper()
    amino_acids = ""

    for i in range(0, len(nucleotide), 3):
        codon = nucleotide[i:i+3]

        if len(codon) != 3:
            amino_acids += "X"
        elif codon == "---":
            amino_acids += "-"
        elif "-" in codon:
            amino_acids += "X"
        else:
            amino_acids += str(Seq(codon).translate())

    protein_sequences.append(amino_acids)

# Calculate conservation
results = []

for position in range(len(protein_sequences[0])):

    column = [seq[position] for seq in protein_sequences]

    # Exclude gaps, unknown amino acids and stop codons
    valid = [aa for aa in column if aa not in ["-", "X", "*"]]

    if not valid:
        continue

    counts = Counter(valid)
    most_common = counts.most_common(1)[0]

    conservation = most_common[1] / len(valid) * 100

    results.append({
        "Position": position + 1,
        "Most_common_amino_acid": most_common[0],
        "Conservation_percent": conservation
    })

# Create output folders
os.makedirs("results/figures", exist_ok=True)
os.makedirs("results/tables", exist_ok=True)

# Save results
df = pd.DataFrame(results)

df.to_csv(
    "results/tables/amino_acid_conservation.csv",
    index=False
)

# Plot conservation
plt.figure(figsize=(14, 5))

plt.plot(
    df["Position"],
    df["Conservation_percent"],
    linewidth=1.5
)

plt.axhline(
    100,
    linestyle="--",
    label="100% conservation"
)

plt.xlabel("HIV-1 Protease Amino Acid Position")
plt.ylabel("Conservation (%)")
plt.title("Amino Acid Conservation Across 50 HIV-1 Protease Sequences")

plt.ylim(0, 105)
plt.xlim(1, 99)

plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()

plt.savefig(
    "results/figures/amino_acid_conservation.png",
    dpi=300
)

plt.show()

print("Conservation analysis completed!")
print("Figure saved successfully.")
print("Table saved successfully.")
