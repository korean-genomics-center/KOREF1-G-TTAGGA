#%%
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
# %%
path_stat_initial = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg.fasta.stats"
path_stat_ntlink = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.ONT_triobin_specific.round5_gapfilling/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg.fa.k32.w100.z1000.ntLink.5rounds.fa.stats"


# %%
table_stat_initial_hap1 = pd.read_csv(path_stat_initial.format(hap = "1"), sep = '\t', names = ["label", "value"])
table_stat_ntlink_hap1 = pd.read_csv(path_stat_ntlink.format(hap = "1"), sep = '\t', names = ["label", "value"])
table_stat_initial_hap2 = pd.read_csv(path_stat_initial.format(hap = "2"), sep = '\t', names = ["label", "value"])
table_stat_ntlink_hap2 = pd.read_csv(path_stat_ntlink.format(hap = "2"), sep = '\t', names = ["label", "value"])
# %%
dict_stat_initial_hap1 = dict(zip(table_stat_initial_hap1["label"], table_stat_initial_hap1["value"]))
dict_stat_ntlink_hap1 = dict(zip(table_stat_ntlink_hap1["label"], table_stat_ntlink_hap1["value"]))
dict_stat_initial_hap2 = dict(zip(table_stat_initial_hap2["label"], table_stat_initial_hap2["value"]))
dict_stat_ntlink_hap2 = dict(zip(table_stat_ntlink_hap2["label"], table_stat_ntlink_hap2["value"]))
# %%
labels_draw = [
    "# scaffolds",
    "# contigs",
    "# gaps",
    "Contig N50"
]

barwidth = 0.2
for label in labels_draw:
    vals_initial = [int(dict_stat_initial_hap1[label]), int(dict_stat_initial_hap2[label])]
    vals_ntlink = [int(dict_stat_ntlink_hap1[label]), int(dict_stat_ntlink_hap2[label])]

    print(1-(vals_ntlink[0] / vals_initial[0]))
    print(1-(vals_ntlink[1] / vals_initial[1]))
    print(1-(sum(vals_ntlink) / sum(vals_initial)))

    plt.figure(figsize=  (2, 2))
    plt.bar([0-barwidth/2, 1-barwidth/2], vals_initial, width = barwidth, color = "gray", label = "Initial Assembly")
    plt.bar([0+barwidth/2, 1+barwidth/2], vals_ntlink, width = barwidth, color = "firebrick", label = "ntLink Scaffolding")
    plt.xticks([0, 1], ["Haplotype 1\nAssembly", "Haplotype 2\nAssembly"])
    plt.title(label)
    # plt.legend(loc = "center left", bbox_to_anchor = (1.01, 0.5))
    plt.show()
    # break
# %%
