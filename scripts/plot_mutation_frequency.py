from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt

reference_record = SeqIO.read(
    "data/sequences/HIV1_protease_99aa.fasta",
    "fasta"
)

reference = str(reference_record.seq).upper()

alignment = list(SeqIO.parse(
    "data/sequences/HIV1_protease_aligned.fasta",
    "fasta"
))

mutation_counts = Counter()

for record in alignment:
    nucleotide_sequence = str(record.seq).upper()
    protein_sequence = str(Seq(nucleotide_sequence).translate())

    if len(protein_sequence) != len(reference):
        continue

    for position, (ref_aa, sample_aa) in enumerate(
        zip(reference, protein_sequence), start=1
    ):
        if ref_aa != sample_aa:
            mutation = f"{ref_aa}{position}{sample_aa}"
            mutation_counts[mutation] += 1

data = []

for mutation, count in mutation_counts.items():
    data.append({
        "Mutation": mutation,
        "Count": count,
        "Frequency (%)": count / len(alignment) * 100
    })

df = pd.DataFrame(data)

# Save complete dataset
df = df.sort_values("Count", ascending=False)
df.to_csv(
    "results/tables/mutation_frequency.csv",
    index=False
)

# Display top 20 substitutions
top20 = df.head(20).sort_values("Frequency (%)")

plt.figure(figsize=(11, 9))

bars = plt.barh(
    top20["Mutation"],
    top20["Frequency (%)"]
)

plt.xlabel("Frequency among 50 sequences (%)")
plt.ylabel("Amino acid substitution")

plt.title(
    "Top 20 Amino Acid Substitutions in HIV-1 Protease"
)

plt.bar_label(
    bars,
    fmt="%.1f%%",
    padding=4
)

plt.xlim(0, top20["Frequency (%)"].max() + 10)

plt.grid(axis="x", linestyle="--", alpha=0.3)
plt.tight_layout()

plt.savefig(
    "results/figures/mutation_frequency.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("Analysis completed!")
print("Total distinct substitutions:", len(df))
print("Displayed substitutions: 20")
print("Figure and complete table saved.")
