from Bio import SeqIO
import random

source = "data/sequences/hiv-db(2).fasta"

records = list(SeqIO.parse(source, "fasta"))

# Keep complete, 297-nt sequences with unambiguous DNA
clean = [
    r for r in records
    if len(r.seq) == 297
    and set(str(r.seq).upper()) <= set("ACGT")
]

# Remove duplicate sequences, retaining the first record
unique = {}
for record in clean:
    sequence = str(record.seq).upper()
    if sequence not in unique:
        unique[sequence] = record

unique_records = list(unique.values())

# Save the full cleaned dataset
SeqIO.write(
    unique_records,
    "data/sequences/HIV1_protease_clean_unique.fasta",
    "fasta"
)

# Reproducible random sample of 50 unique sequences
rng = random.Random(42)
sample = rng.sample(unique_records, 50)

SeqIO.write(
    sample,
    "data/sequences/HIV1_protease_sample50.fasta",
    "fasta"
)

print("Total downloaded:", len(records))
print("Clean sequences:", len(clean))
print("Unique sequences:", len(unique_records))
print("Selected for analysis:", len(sample))
print("Random seed: 42")

