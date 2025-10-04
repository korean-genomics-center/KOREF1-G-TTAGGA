#%%
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

#%%
path_hamming = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20.with_0/Hap1/Merqury_Hamming/assembly.hapmers.count"

col_origin = "PATERNAL_Trimmed_F15_T5.noNbase.merged.k21.hapmer"
#%%
table_hamming = pd.read_csv(path_hamming, sep = '\t')
# %%
table_hamming["tot_ham"] = table_hamming["PATERNAL_Trimmed_F15_T5.noNbase.merged.k21.hapmer"] + table_hamming["MATERNAL_Trimmed_F15_T5.noNbase.merged.k21.hapmer"]

table_hamming["hamming_error"] = 1 - table_hamming[col_origin]/table_hamming["tot_ham"]
# %%
plt.rcParams["font.size"] = 15
plt.figure(figsize = (5, 4))
plt.hist(table_hamming["hamming_error"], bins = np.linspace(0, 1, 101), color = "gray")
plt.xlabel("Hamming Error")
plt.ylabel("# Contigs")
plt.show()
# %%
plt.rcParams["font.size"] = 15
plt.figure(figsize = (5, 4))
plt.scatter(table_hamming["Size"], table_hamming["tot_ham"], color = "gray", s = 3)
plt.xlabel("Contig size (bp)")
plt.ylabel("# Total Hapmers")
# plt.ylim(bottom = -1e4, top = 1e5)
plt.show()
# %%
plt.rcParams["font.size"] = 15
plt.figure(figsize = (5, 4))
plt.scatter(table_hamming["Size"], table_hamming["hamming_error"], color = "gray", s = 3)
plt.xlabel("Contig size (bp)")
plt.ylabel("Hamming Error")
# plt.xlim(left = -1e4, right = 1e6)
plt.show()

# %%
plt.rcParams["font.size"] = 15
plt.figure(figsize = (5, 4))
plt.scatter(table_hamming["tot_ham"], table_hamming["hamming_error"], color = "gray", s = 3)
plt.xlabel("# Total Hapmers")
plt.ylabel("Hamming Error")
plt.xlim(left = -1e3, right = 1e4)
plt.show()

# %%
print("Hamming Error < 0.05: ", sum(table_hamming["hamming_error"] < 0.05) / table_hamming.shape[0])
# %%
print("Hapmer count == 0: ",sum(table_hamming["tot_ham"] == 0) / table_hamming.shape[0])
# %%
print("Total Hamming Error: ", (1 - table_hamming[col_origin].sum()/table_hamming["tot_ham"].sum()) * 100, "%")
# %%
