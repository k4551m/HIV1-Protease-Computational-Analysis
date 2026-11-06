import pandas as pd
import matplotlib.pyplot as plt

# Load Stanford HIVdb results
df = pd.read_csv("data/resistanceSummaries.csv")

# Protease inhibitors and their corresponding columns
drugs = {
    "Atazanavir/r": "ATV/r Level",
    "Darunavir/r": "DRV/r Level",
    "Fosamprenavir/r": "FPV/r Level",
    "Indinavir/r": "IDV/r Level",
    "Lopinavir/r": "LPV/r Level",
    "Nelfinavir": "NFV Level",
    "Saquinavir/r": "SQV/r Level",
    "Tipranavir/r": "TPV/r Level"
}

# HIVdb interpretation categories
categories = {
    1: "Susceptible",
    2: "Potential low-level",
    3: "Low-level",
    4: "Intermediate",
    5: "High-level"
}

# Calculate category counts
results = {}

for drug, column in drugs.items():
    counts = (
        df[column]
        .value_counts()
        .reindex(categories.keys(), fill_value=0)
    )

    results[drug] = counts.values

# Create results table
summary = pd.DataFrame(
    results,
    index=categories.values()
).T

summary.to_csv(
    "results/tables/drug_resistance_summary.csv"
)

# Create stacked horizontal bar chart
ax = summary.plot(
    kind="barh",
    stacked=True,
    figsize=(12, 7),
    width=0.7
)

# Add count labels to visible segments
for container in ax.containers:
    labels = [
        str(int(value)) if value >= 3 else ""
        for value in container.datavalues
    ]

    ax.bar_label(
        container,
        labels=labels,
        label_type="center",
        fontsize=8
    )

plt.xlabel("Number of sequences")
plt.ylabel("Protease inhibitor")

plt.title(
    "HIV-1 Protease Inhibitor Resistance Predictions"
)

plt.legend(
    title="Resistance interpretation",
    bbox_to_anchor=(1.02, 1),
    loc="upper left"
)

plt.xlim(0, len(df))
plt.tight_layout()

plt.savefig(
    "results/figures/drug_resistance_profile.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("Drug resistance analysis completed!")
print("Sequences analysed:", len(df))
print("Protease inhibitors analysed:", len(drugs))
print("\nResistance summary:")
print(summary)
print("\nFigure and table saved successfully.")
