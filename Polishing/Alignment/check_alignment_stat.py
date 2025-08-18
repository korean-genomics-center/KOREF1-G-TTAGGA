#%%
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
# %%
table_identity_hap1_p = pd.read_csv("/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Read_Alignment/Align.Winnowmap2/Align_Gapclosed_ASM/Primary_Alignment/Winnowmap_Align.Primary.Ref.KOREF_ASM.hap1.QRY.PacBio_HiFi_ccs_converted.Filter_L1000_QV20.triobin_split.triobin_p.mapq_identity", sep = ' ', names = ["readname", "mapq", "alignlen", "identity"])
table_identity_hap1_m = pd.read_csv("/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Read_Alignment/Align.Winnowmap2/Align_Gapclosed_ASM/Primary_Alignment/Winnowmap_Align.Primary.Ref.KOREF_ASM.hap1.QRY.PacBio_HiFi_ccs_converted.Filter_L1000_QV20.triobin_split.triobin_m.mapq_identity", sep = ' ', names = ["readname", "mapq", "alignlen", "identity"])
# %%
plt.rcParams["font.size"] = 15
plt.hist(
    table_identity_hap1_p["mapq"],
    bins = np.linspace(0, 60, 61),
    color = "gray"
)
# plt.yscale("log")
plt.xlabel("Mapping Quality")
plt.ylabel("Alignments")
plt.show()
# %%
plt.rcParams["font.size"] = 15
plt.hist(
    table_identity_hap1_m["mapq"],
    bins = np.linspace(0, 60, 61),
    color = "gray"
)
# plt.yscale("log")
plt.xlabel("Mapping Quality")
plt.ylabel("Alignments")
plt.show()

# %%
plt.rcParams["font.size"] = 15
plt.hist(
    table_identity_hap1_p["identity"],
    bins = np.linspace(0, 100, 101),
    color = "gray"
)
# plt.yscale("log")
plt.xlabel("Identity (%)")
plt.ylabel("Alignments")
plt.show()

# %%
plt.rcParams["font.size"] = 15
plt.hist(
    table_identity_hap1_m["identity"],
    bins = np.linspace(0, 100, 101),
    color = "gray"
)
# plt.yscale("log")
plt.xlabel("Identity (%)")
plt.ylabel("Alignments")
plt.show()

# %%
