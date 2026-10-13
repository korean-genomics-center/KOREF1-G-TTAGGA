#%%
import os, math

import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np

#%%
path_bedpe = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Repeat_Annotation/Segmental_Duplication/Biser/Hap2/FinalASM.KOREF_S1.Trio_Specific.hap2.Reordered.masked.SD.check_overlap"

thresh_filt = 50
path_save = f"{path_bedpe}.filt_overlap_{thresh_filt}"

#%%
cols_original = ["chr1", "start1", "end1", "chr2", "start2", "end2", "reference", "score", "strand1", "strand2", "max_len", "aln_len", "cigar", "optional"]
#%%
table_bedpe = pd.read_csv(
    path_bedpe, sep = '\t', names = cols_original + ["overlap_len"]
)

# %%
table_bedpe["region1_len"] = table_bedpe["end1"] - table_bedpe["start1"]
table_bedpe["region2_len"] = table_bedpe["end2"] - table_bedpe["start2"]

table_bedpe["overlap_ratio"] = table_bedpe["overlap_len"] / table_bedpe[["region1_len", "region2_len"]].max(axis = 1)
# %%
plt.rcParams["font.size"] = 14
plt.rcParams["font.family"] = "sans-serif"

plt.figure(figsize = (4, 3))
plt.hist(
    table_bedpe["overlap_ratio"], bins = np.linspace(0, 1, 101), color = "gray"
)
plt.xlabel("Ratio of Overlap")
plt.ylabel("Count")
# plt.yscale("log")
plt.show()
# %%
table_bedpe_filt = table_bedpe[table_bedpe["overlap_ratio"] * 100 < thresh_filt]

#%%
def get_bedpe_pair_id(chr1, start1, end1, chr2, start2, end2):
    if chr1 == chr2:
        if start1 < start2:
            return f"{chr1}:{start1}:{end1}__{chr2}:{start2}:{end2}"
        else:
            return f"{chr2}:{start2}:{end2}__{chr1}:{start1}:{end1}"
    else:
        if chr1 < chr2:
            return f"{chr1}:{start1}:{end1}__{chr2}:{start2}:{end2}"
        else:
            return f"{chr2}:{start2}:{end2}__{chr1}:{start1}:{end1}"

table_bedpe_filt["relid"] = table_bedpe_filt.apply(lambda row: get_bedpe_pair_id(row["chr1"], row["start1"], row["end1"], row["chr2"], row["start2"], row["end2"]), axis = 1)
#%%
table_bedpe_filt_unique = table_bedpe_filt.sort_values(by = ["score"]).drop_duplicates(subset = ["relid"])
# %%
table_bedpe_filt_unique.sort_values(by = ["chr1", "start1"])[cols_original].to_csv(path_save, sep = '\t', header = False, index = False)
# %%
