#%%
import math, sys, os
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

sys.path.append(str(Path(__file__).parents[3] / "common"))
from fasta_parser import read_fasta_as_dict

#%%
hap = 2
triobin = 'm'

path_fasta = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/ONT_UL.Filter_L50kb_QV20_Q30Perc80/Hap{hap}/assembly.fasta"
# path_qv = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.Binned_Input/ONT_UL.Filter_L50kb_QV20_Q30Perc80.HiFi.Filter_L1000_QV20/Hap{hap}/Merqury_QV/KOREF_S1.Binned_Assembly.ul_L50kb_QV20_Q30Perc80.triobin_{triobin}.hifi_L1000_QV20.triobin_{triobin}.haploid.l0.asm.bp.p_ctg.qv"
path_qv = os.path.join(
    os.path.dirname(path_fasta),
    "Merqury_QV",
    f"{os.path.basename(path_fasta).replace('.fasta', '')}.perseq.qv"
)

#%%
table_qv = pd.read_csv(path_qv, names = ["seqname", "wrong", "total", "qv", "error"], sep = '\t')

dict_fasta = read_fasta_as_dict(path_fasta)

dict_seqname_to_len = dict(zip(
    table_qv["seqname"].to_list(),
    list(map(lambda name: len(dict_fasta[name]), table_qv["seqname"].to_list()))
))
# %%
table_qv["qv_replace"] = table_qv["qv"].apply(lambda val:99 if np.isinf(val) else val)
# %%
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.size"] = 15
#%%
plt.figure(figsize = (5, 4))
plt.hist(table_qv["qv_replace"], bins = np.linspace(0, 100, 101), color = "gray")
plt.xlabel("QV per Contig")
plt.ylabel("# Contigs")
plt.show()
# %%
plt.figure(figsize = (5, 4))
plt.scatter(
    list(map(dict_seqname_to_len.__getitem__, table_qv["seqname"])),
    table_qv["qv_replace"],
    color = "gray",
    s = 3
)
plt.xlabel("Contig Length (bp)")
plt.ylabel("Contig QV")
plt.show()

# %%
