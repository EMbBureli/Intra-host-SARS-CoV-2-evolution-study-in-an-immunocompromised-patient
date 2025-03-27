# génération du fichier résultat : pN_pN_mutations_final.csv, qui comporte les colonnes : Echantillon	Region	Mutations_Synonymes_iSNV_WH	Mutations_Non_Synonymes_iSNV_WH	Mutations_Synonymes_iSNV_TOT	Mutations_Non_Synonymes_iSNV_TOT	Mutations_Synonymes_SNP_WH	Mutations_Non_Synonymes_SNP_WH	Mutations_Synonymes_SNP_TOT	Mutations_Non_Synonymes_SNP_TOT

import pandas as pd
import numpy as np
import os

MUTATION_FILES = {
    "iSNV_WH": "filtered_iSNV_WH_50x.csv",
    "iSNV_TOT": "filtered_iSNV_TOT_50x.csv",
    "SNP_WH": "filtered_SNP_WH_50x.csv",
    "SNP_TOT": "filtered_SNP_TOT_50x.csv"
}

GENOMIC_REGIONS = [
    ("nsp1", 266, 805), ("nsp2", 806, 2719), ("nsp3", 2720, 8554), ("nsp4", 8555, 10054),
    ("nsp5", 10055, 10972), ("nsp6", 10973, 11842), ("nsp7", 11843, 12091), ("nsp8", 12092, 12685),
    ("nsp9", 12686, 13024), ("nsp10", 13025, 13441), ("nsp11", 13442, 13468), ("nsp12 (RdRp)", 13468, 16236),
    ("nsp13", 16237, 18039), ("nsp14", 18040, 19620), ("nsp15", 19621, 20658), ("nsp16", 20659, 21555),
    ("S", 21563, 25384), ("ORF3a", 25393, 26220), ("E", 26245, 26472), ("M", 26523, 27191),
    ("ORF6", 27202, 27387), ("ORF7a", 27394, 27759), ("ORF7b", 27756, 27887), ("ORF8", 27894, 28259),
    ("N", 28274, 29533), ("ORF9b", 28284, 28577), ("ORF10", 29558, 29674)
]

#  séparer les chevauchements
def assign_region(pos):
    """
    Renvoie une liste de régions si une mutation appartient à plusieurs ORFs.
    Exemple :
        - POS = 28284 → ["N", "ORF9b"] (car chevauchement)
        - POS = 28000 → ["N"] (pas de chevauchement)
    """
    regions = [name for name, start, end in GENOMIC_REGIONS if start <= pos <= end]
    return regions if regions else ["Unknown"]

# classifie les mutations synonymes et non-synonymes
def classify_mutations(df):
    """
    Ajoute deux colonnes :
    - "Mutation_Synonyme" (1 si synonyme, 0 sinon)
    - "Mutation_Non_Synonyme" (1 si non-synonyme, 0 sinon)
    Exclut les mutations non codantes (NA).
    """
    df["Mutation_Synonyme"] = np.where((df["ALT_AA"] == df["REF_AA"]) & (df["ALT_AA"].notna()), 1, 0)
    df["Mutation_Non_Synonyme"] = np.where((df["ALT_AA"] != df["REF_AA"]) & (df["ALT_AA"].notna()), 1, 0)
    return df
results = []


for key, file in MUTATION_FILES.items():
    if os.path.exists(file):
        df = pd.read_csv(file, sep=";", encoding="utf-8")
        df.columns = df.columns.str.strip()  # Supprime les espaces cachés

        df["Regions"] = df["POS"].apply(assign_region)
        df_exploded = df.explode("Regions").rename(columns={"Regions": "Region"})

        # Classification des mutations synonymes et non-synonymes
        df_exploded = classify_mutations(df_exploded)

        # Agrégation par échantillon et région
        df_grouped = df_exploded.groupby(["Echantillon", "Region"])[["Mutation_Synonyme", "Mutation_Non_Synonyme"]].sum().reset_index()

    
        df_grouped.rename(columns={
            "Mutation_Synonyme": f"Mutations_Synonymes_{key}",
            "Mutation_Non_Synonyme": f"Mutations_Non_Synonymes_{key}"
        }, inplace=True)

        results.append(df_grouped)
    else:
        print(f"il manque le fichier: {file}")

# Fusionne les résultats pour les 4 classes de mutations
df_final = results[0]
for df in results[1:]:
    df_final = df_final.merge(df, on=["Echantillon", "Region"], how="outer")

#  les NaN sont remplacés par 0 et convertis en int
df_final.fillna(0, inplace=True)
df_final[df_final.select_dtypes(include=[np.number]).columns] = df_final.select_dtypes(include=[np.number]).astype(int)

output_file = "pN_pS_mutations_final.csv"
df_final.to_csv(output_file, sep=";", encoding="utf-8", index=False)
