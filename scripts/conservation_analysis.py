from Bio import AlignIO
from collections import Counter

alignment = AlignIO.read(
    "data/sequences/HIV1_protease_aligned.fasta",
    "fasta"
)

print("Sequences:", len(alignment))
print("Alignment length:", alignment.get_alignment_length())

for i in range(alignment.get_alignment_length()):
    column = alignment[:, i]
    counts = Counter(column)
    most_common = counts.most_common(1)[0]

    conservation = (most_common[1] / len(alignment)) * 100

    print(
        f"Position {i+1}: "
        f"Most common = {most_common[0]}, "
        f"Conservation = {conservation:.1f}%"
    )
