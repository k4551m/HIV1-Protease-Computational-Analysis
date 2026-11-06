from Bio import AlignIO
from collections import Counter

alignment = AlignIO.read(
    "data/sequences/HIV1_protease_aligned.fasta",
    "fasta"
)

protein_sequences = []

for record in alignment:
    protein = str(record.seq).upper().replace("U", "T")

    # Translate nucleotide sequence
    amino_acids = ""

    for i in range(0, len(protein), 3):
        codon = protein[i:i+3]
        from Bio.Seq import Seq
        amino_acids += str(Seq(codon).translate())

    protein_sequences.append(amino_acids)

print("Number of protein sequences:", len(protein_sequences))
print("Protein length:", len(protein_sequences[0]))

print("\nVARIABLE AMINO-ACID POSITIONS")

for i in range(len(protein_sequences[0])):
    column = [seq[i] for seq in protein_sequences]
    counts = Counter(column)

    most_common = counts.most_common(1)[0]
    conservation = most_common[1] / len(column) * 100

    if conservation < 100:
        print(
            f"Position {i+1}: "
            f"Most common = {most_common[0]}, "
            f"Conservation = {conservation:.1f}%"
        )

print("\nQUALITY CONTROL")

internal_stops = sum("*" in protein[:-1] for protein in protein_sequences)

print("Sequences with internal stop codons:", internal_stops)
print("Expected protein length: 99 amino acids")
print("Actual protein length:", len(protein_sequences[0]))
