import pandas as pd
import os

INPUT_FILES = [
    "D0.csv", "D229.csv", "D276.csv", "D277.csv", "D285.csv", "D290.csv",
    "D329.csv", "D330.csv", "D378.csv", "D431.csv", "D432.csv", "D434.csv"
]
MARKERS_FILE = "mut_markers.csv"  #mutations marqueurs à exclure

#délétions à exclure
DELETION_POSITIONS = {21990, 21764}

def load_markers(file):
    if os.path.exists(file):
        markers_df = pd.read_csv(file, sep=";")
        if {"POS", "REF", "ALT"}.issubset(markers_df.columns):
            return set(zip(markers_df["POS"], markers_df["REF"], markers_df["ALT"]))
    return set()

EXCLUDED_MARKERS = load_markers(MARKERS_FILE)

# Dictionnaire pour stocker les résultats de chaque filtre
results = {
    "iSNV_WH": [],
    "iSNV_TOT": [],
    "SNP_WH": [],
    "SNP_TOT": []
}

#boucle pour filtrer
for file in INPUT_FILES:
    if not os.path.exists(file):
        print(f"Fichier introuvable : {file}, ignoré.")
        continue

    df = pd.read_csv(file, sep=";", engine="python")

    #Colonnes requises
    required_columns = {"POS", "REF", "ALT", "TOTAL_DP", "ALT_FREQ", "REF_CODON", "ALT_CODON", "REF_AA", "ALT_AA"}
    if not required_columns.issubset(df.columns):
        print(f" Fichier {file} ignoré - manque colonnes : {required_columns - set(df.columns)}")
        continue

    df_filtered = df[df["TOTAL_DP"] >= 50]

    # Filtrage pour pool des iSNV (FA >= 0.05 et < 0.5)
    df_iSNV = df_filtered[(df_filtered["ALT_FREQ"] >= 0.05) & (df_filtered["ALT_FREQ"] < 0.5)]
    df_iSNV_WH = df_iSNV[~df_iSNV.apply(lambda row: (row["POS"], row["REF"], row["ALT"]) in EXCLUDED_MARKERS, axis=1)]

    #Filtrage pour SNP (Fréquence allélique >= 0.5)
    df_SNP = df_filtered[df_filtered["ALT_FREQ"] >= 0.5]
    df_SNP_WH = df_SNP[~df_SNP.apply(lambda row: (row["POS"], row["REF"], row["ALT"]) in EXCLUDED_MARKERS, axis=1)]

    # Ajoute le nom du fichier.csv comme échantillon
    for df_sub, key in zip([df_iSNV_WH, df_iSNV, df_SNP_WH, df_SNP], ["iSNV_WH", "iSNV_TOT", "SNP_WH", "SNP_TOT"]):
        df_sub["Echantillon"] = file.replace(".csv", "")
        results[key].append(df_sub[["Echantillon", "POS", "REF", "ALT", "ALT_FREQ", "REF_CODON", "REF_AA", "ALT_CODON", "ALT_AA"]])

output_files = {
    "iSNV_WH": "filtered_iSNV_WH.csv",
    "iSNV_TOT": "filtered_iSNV_TOT.csv",
    "SNP_WH": "filtered_SNP_WH.csv",
    "SNP_TOT": "filtered_SNP_TOT.csv"
}

for key, file_name in output_files.items():
    if results[key]:
        df_results = pd.concat(results[key], ignore_index=True)
        df_results.to_csv(file_name, sep=";", index=False)
        print(f"\n Done, Résultats exportés dans '{file_name}'.")
    else:
        print(f"\n No result for '{file_name}'")
