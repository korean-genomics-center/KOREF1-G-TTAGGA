#%%
import os, sys
from pathlib import Path

import pandas as pd

sys.path.append(str(Path(__file__).parents[2]))
from fasta_parser import *
# %%
path_fasta1 = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Final_Assembly/FinalASM.KOREF_S1.Trio_Specific.hap1.Reordered.fasta"
path_fasta2 = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Final_Assembly/FinalASM.KOREF_S1.Trio_Specific.hap2.Reordered.fasta"

path_vote = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Final_Assembly/Vote_Chromosomes_for_Combined.tsv"

path_save_combined = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Final_Assembly/FinalASM.KOREF_S1.Trio_Specific.Combined.Reordered.fasta"
# %%
dict_fasta1 = read_fasta_as_dict(path_fasta1)
dict_fasta2 = read_fasta_as_dict(path_fasta2)
table_vote = pd.read_csv(path_vote, sep = '\t', index_col = 0)
# %%
dict_combined_fasta = dict()
for chrname, selected in table_vote["Combined"].to_dict().items():
    if selected == "Paternal":
        dict_combined_fasta[chrname] = dict_fasta1[chrname]
    elif selected == "Maternal":
        dict_combined_fasta[chrname] = dict_fasta2[chrname]
    else:
        raise Exception(f"Unexpected Chromosome:", chrname)
# %%
save_dict_fasta_into_file(dict_combined_fasta, path_save_combined)
# %%
