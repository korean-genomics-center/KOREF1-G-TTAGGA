#%%
import math

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
# %%
hap = 1

path_qv_total_result = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Reference_Fasta/Round2.SNV_Polished_Fasta/Merqury/PolishedASM.KOREF_S1.Polish_R2.SNV.merqury.k21/PolishedASM.KOREF_S1.Polish_R2.SNV.merqury.k21.PolishedASM.KOREF_S1.Trio_Specific.UpperCase.Polish_R2.SNV.hap{hap}.qv"
path_qv_per_window_result = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Reference_Fasta/Round2.SNV_Polished_Fasta/Merqury/PolishedASM.KOREF_S1.Polish_R2.SNV.QV_per_window.k21/PolishedASM.KOREF_S1.Polish_R2.SNV.Hap{hap}.QV_per_window.k21.chunked.100000.qv"

#%%
table_qv_total = pd.read_csv(path_qv_total_result, names = ["seqname", "wrong", "total", "qv", "error"], sep = '\t')
table_qv = pd.read_csv(path_qv_per_window_result, names = ["seqname", "wrong", "total", "qv", "error"], sep = '\t')
# %%
table_qv.loc[np.isinf(table_qv["qv"]), "qv"] = 99
# %%
table_qv["chrname"] = table_qv["seqname"].apply(lambda val: val.split(':')[-2])
table_qv["range_start"] = table_qv["seqname"].apply(lambda val: int(val.split(':')[-1].split('-')[0]))
table_qv["range_end"] = table_qv["seqname"].apply(lambda val: int(val.split(':')[-1].split('-')[1]))
table_qv["range_mid"] = (table_qv["range_start"]+table_qv["range_end"]) / 2
# %%
list_chrname = list(map(lambda val: f"chr{val}", range(1, 23))) + ["chrX", "chrY"]
list_chrname = list(filter(lambda val: val in table_qv["chrname"].to_list(), list_chrname))
# %%
dict_chrname_to_maxlen = dict()
for chrname in list_chrname:
    dict_chrname_to_maxlen[chrname] = table_qv[table_qv["chrname"] == chrname]["range_end"].max()
dict_chrname_to_qvtot = dict(zip(table_qv_total["seqname"], table_qv_total["qv"]))
# %%
plt.rcParams['font.family'] = 'Sans-serif'
plt.rcParams['font.size'] = 15

ax_height = 100
ax_height_pad = 50

total_qv_plot_ratio = 0.4
ax_width_pad_ratio = 0.075

length_divider = 1e4
xlim_padding = 1e6
qv_alert = 40

col_length_of_scatter = math.ceil(max(dict_chrname_to_maxlen.values())/length_divider)

fig = plt.figure(figsize=(12, 15))
row = (ax_height+ax_height_pad) * len(list_chrname) + ax_height
col = math.ceil(col_length_of_scatter*(1+total_qv_plot_ratio+ax_width_pad_ratio))

totalqv_width = math.floor(col_length_of_scatter*total_qv_plot_ratio)
ax_width_pad = math.floor(col_length_of_scatter*ax_width_pad_ratio)

# fontsize_label = 17
gsfig = GridSpec(
    row, col, 
    left=0, right=1, bottom=0,
    top=1, wspace=1, hspace=1)

dict_ax = dict()
for ind, chrname in enumerate(list_chrname):
    gsrange = gsfig[(ax_height+ax_height_pad)*ind:(ax_height+ax_height_pad)*ind+ax_height, totalqv_width+ax_width_pad:totalqv_width+ax_width_pad+math.ceil(dict_chrname_to_maxlen[chrname]/length_divider)]
    dict_ax[chrname] = fig.add_subplot(gsrange)

for chrname in list_chrname:
    table_qv_chr = table_qv[table_qv["chrname"] == chrname]

    color = list(map(lambda val: "gray" if val > qv_alert else "r", table_qv_chr["qv"]))
    size = list(map(lambda val: 10 if val > qv_alert else 3, table_qv_chr["qv"]))
    dict_ax[chrname].scatter(
        table_qv_chr["range_mid"], table_qv_chr["qv"],
        s = size, color = color, zorder = 5
    )
    dict_ax[chrname].scatter(
        table_qv_chr[table_qv_chr["qv"]<=qv_alert]["range_mid"], table_qv_chr[table_qv_chr["qv"]<=qv_alert]["qv"],
        s = 100, color = 'r', zorder = 7, marker = '*'
    )
    dict_ax[chrname].plot(
        table_qv_chr["range_mid"], table_qv_chr["qv"],
        linewidth = 0.75, color = "gray", zorder = 3
    )

    dict_ax[chrname].set_ylim(0, 110)
    dict_ax[chrname].set_xlim(0-xlim_padding, dict_chrname_to_maxlen[chrname]+xlim_padding)
    dict_ax[chrname].set_xticks([])
    dict_ax[chrname].axhline(qv_alert, linewidth = 1, linestyle = "--", color = "gray", zorder = 1)
    # dict_ax[chrname].set_ylabel("QV")
    dict_ax[chrname].annotate(chrname, (1+0.01 * max(dict_chrname_to_maxlen.values())/dict_chrname_to_maxlen[chrname], 0.5), va = "center", ha = "left", xycoords = "axes fraction", fontsize = plt.rcParams["font.size"]+5, fontweight = "bold")

gsrange = gsfig[0:(ax_height+ax_height_pad)*(len(list_chrname)-1)+ax_height, 0:totalqv_width]
ax_bar = fig.add_subplot(gsrange)
ax_bar.barh(
    list_chrname[::-1], list(map(lambda name: 99 if np.isinf(dict_chrname_to_qvtot[name]) else dict_chrname_to_qvtot[name], list_chrname))[::-1],
    color = "gray"
)
for ind, chrname in enumerate(list_chrname[::-1]):
    ax_bar.annotate(f"{dict_chrname_to_qvtot[chrname]:.2f}", (dict_chrname_to_qvtot[chrname], ind), va = "center", ha = "right", fontsize = plt.rcParams["font.size"]+3)
ax_bar.invert_xaxis()
ax_bar.set_yticks([])
ax_bar.set_ylim([-0.5, len(list_chrname)-0.5])
ax_bar.set_xlim(110, 0)
plt.show()
# %%
fig, ax = plt.subplots(1, 1, figsize = (7, 3))
chrname = "chr22"
table_qv_chr = table_qv[table_qv["chrname"] == chrname]

color = list(map(lambda val: "gray" if val > qv_alert else "firebrick", table_qv_chr["qv"]))
ax.scatter(
    table_qv_chr["range_mid"], table_qv_chr["qv"],
    s = 5, color = color, zorder = 5
)
ax.plot(
    table_qv_chr["range_mid"], table_qv_chr["qv"],
    linewidth = 1, color = "gray", zorder = 3
)

ax.set_ylim(0, 110)
ax.set_xlim(0-xlim_padding, dict_chrname_to_maxlen[chrname]+xlim_padding)
# ax.set_xticks([])
plt.xlabel("Coordinate (bp)")
plt.ylabel("QV")
ax.axhline(qv_alert, linewidth = 1, linestyle = "--", color = "gray", zorder = 1)
plt.show()
# %%
