from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter

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

    protein_sequence = str(
        Seq(nucleotide_sequence).translate()
    )

    if len(protein_sequence) != len(reference):
        print("Length mismatch:", record.id)
        continue

    for position, (ref_aa, sample_aa) in enumerate(
        zip(reference, protein_sequence), start=1
    ):
        if ref_aa != sample_aa:
            mutation = f"{ref_aa}{position}{sample_aa}"
            mutation_counts[mutation] += 1

print("REFERENCE LENGTH:", len(reference))
print("SEQUENCES ANALYSED:", len(alignment))

print("\nIDENTIFIED AMINO-ACID SUBSTITUTIONS")

for mutation, count in sorted(
    mutation_counts.items(),
    key=lambda item: (int(item[0][1:-1]), item[0])
):
    frequency = count / len(alignment) * 100

    print(
        f"{mutation}: {count}/{len(alignment)} "
        f"sequences ({frequency:.1f}%)"
    )

print("\nTOTAL DISTINCT SUBSTITUTIONS:", len(mutation_counts))
