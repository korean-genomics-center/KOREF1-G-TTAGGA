#%%
import os, math

import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np

#%%
path_bedpe = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Repeat_Annotation/Segmental_Duplication/Biser/Hap1/FinalASM.KOREF_S1.Trio_Specific.hap1.Reordered.masked.SD.check_overlap.filt_overlap_50"

thresh_err = 10
path_save = f"{path_bedpe}.err_lt{thresh_err}"

#%%
cols_original = ["chr1", "start1", "end1", "chr2", "start2", "end2", "reference", "score", "strand1", "strand2", "max_len", "aln_len", "cigar", "optional"]
# %%
table_bedpe = pd.read_csv(
    path_bedpe, sep = '\t', names = cols_original
)
# %%
plt.rcParams["font.size"] = 14
plt.rcParams["font.family"] = "sans-serif"

plt.figure(figsize = (4, 3))
plt.hist(
    table_bedpe["score"], bins = np.linspace(0, 100, 101), color = "gray"
)
plt.xlabel("Alignment Error (%)")
plt.ylabel("Count")
# plt.yscale("log")
plt.axvline(thresh_err, linewidth = 2, color = "firebrick", linestyle = "--")
plt.show()
# %%
plt.rcParams["font.size"] = 14
plt.rcParams["font.family"] = "sans-serif"

plt.figure(figsize = (4, 3))
plt.hist(
    table_bedpe["max_len"], bins = 100, color = "gray"
)
plt.xlabel("Length of Mate (bp)")
plt.ylabel("Count")
# plt.yscale("log")
# plt.axvline(thresh_err, linewidth = 2, color = "firebrick", linestyle = "--")
plt.show()
# %%
plt.rcParams["font.size"] = 14
plt.rcParams["font.family"] = "sans-serif"

plt.figure(figsize = (4, 3))
plt.hist(
    table_bedpe["max_len"], bins = 100, color = "gray"
)
plt.xlabel("Length of Mate (bp)")
plt.ylabel("Count")
plt.yscale("log")
# plt.axvline(thresh_err, linewidth = 2, color = "firebrick", linestyle = "--")
plt.show()
# %%
table_bedpe_filt = table_bedpe[table_bedpe["score"] < thresh_err]
# %%
table_bedpe_filt[cols_original].to_csv(path_save, sep = '\t', header = False, index = False)
# %%
