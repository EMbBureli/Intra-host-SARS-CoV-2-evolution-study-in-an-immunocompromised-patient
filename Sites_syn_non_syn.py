from Bio import SeqIO
from Bio.Data import CodonTable
import pandas as pd
import os

REFERENCE_FASTA = "Wuhan.fasta"

# Définition des régions du génome en ignorant le codon start et le codon stop pour ne prendre en compte ques les séquences codantes
GENOME_REGIONS = {
    "orf1a": (266 + 3, 13468 - 3), "nsp1": (266 + 3, 805 - 3), "nsp2": (806 + 3, 2719 - 3),
    "nsp3": (2720 + 3, 8554 - 3), "nsp4": (8555 + 3, 10054 - 3), "nsp5": (10055 + 3, 10972 - 3),
    "nsp6": (10973 + 3, 11842 - 3), "nsp7": (11843 + 3, 12091 - 3), "nsp8": (12092 + 3, 12685 - 3),
    "nsp9": (12686 + 3, 13024 - 3), "nsp10": (13025 + 3, 13441 - 3), "nsp11": (13442 + 3, 13468 - 3),
    "orf1b": (13468 + 3, 21555 - 3), "nsp12": (13468 + 3, 16236 - 3), "nsp13": (16237 + 3, 18039 - 3),
    "nsp14": (18040 + 3, 19620 - 3), "nsp15": (19621 + 3, 20658 - 3), "nsp16": (20659 + 3, 21555 - 3),
    "S": (21563 + 3, 25384 - 3), "ORF3a": (25393 + 3, 26220 - 3), "E": (26245 + 3, 26472 - 3),
    "M": (26523 + 3, 27191 - 3), "ORF6": (27202 + 3, 27387 - 3), "ORF7a": (27394 + 3, 27759 - 3),
    "ORF7b": (27756 + 3, 27887 - 3), "ORF8": (27894 + 3, 28259 - 3), "N": (28274 + 3, 29533 - 3),
    "ORF9b": (28284 + 3, 28577 - 3), "ORF10": (29558 + 3, 29674 - 3),
}


def load_reference_sequence(fasta_file):
    if not os.path.exists(fasta_file):
        print(f" Erreur : Le fichier {fasta_file} n'a pas été trouvé")
        exit()

    record = next(SeqIO.parse(fasta_file, "fasta"))
    return str(record.seq).upper()

reference_seq = load_reference_sequence(REFERENCE_FASTA)

# code génétique
standard_table = CodonTable.unambiguous_dna_by_id[1].forward_table

# Calcul des sites synonymes et non-synonymes
def compute_syn_nonsyn_sites():
    site_info = []

    for region, (start, end) in GENOME_REGIONS.items():
        syn_sites, nonsyn_sites = 0, 0

        # Vérification pour éviter des erreurs de séquence
        if start >= end:
            print(f"Région {region} invalide (start >= end), ignorée.")
            continue

        # Analyse des sites synonymes et non-synonymes
        for pos in range(start, end, 3):
            codon = reference_seq[pos - 1:pos + 2]  # Extraction du codon

            if len(codon) == 3 and codon in standard_table:
                aa = standard_table[codon]  # Acide aminé correspondant dans la table

                for i in range(3):
                    for base in "ATGC":
                        if base != codon[i]:  # Mutation sur une seule base
                            mutated_codon = codon[:i] + base + codon[i + 1:]

                            if mutated_codon in standard_table:
                                if standard_table[mutated_codon] == aa:
                                    syn_sites += 1 / 3  # Mutation synonyme
                                else:
                                    nonsyn_sites += 1 / 3  # Mutation non-synonyme

        site_info.append({"Région": region, "Sites_Synonymes": syn_sites, "Sites_Non_Synonymes": nonsyn_sites})

    return pd.DataFrame(site_info)

# Calcul des sites pour toutes les régions 
df_sites = compute_syn_nonsyn_sites()

output_filename = "sites_syn_nonsyn.csv"
df_sites.to_csv(output_filename, sep=";", index=False)
print(df_sites)
