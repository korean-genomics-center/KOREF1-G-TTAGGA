#%%
import math

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
# %%
k=21
window = 10000
#%%
path_qv_total_result_format = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Merqury_Evaluation/Merqury_QV/PolishedASM.KOREF_S1.Trio_Specific.UpperCase.Polish_R8.SV_SNV.hap{{hap}}.k{k}.perseq.qv"
path_qv_per_window_result_format = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Merqury_Evaluation/Merqury_QV_per_window/PolishedASM.KOREF_S1.Trio_Specific.UpperCase.Polish_R8.SV_SNV.hap{{hap}}.QV_per_window.k{k}.chunked.{window}.qv"

#%%
dict_table_qv_total = dict()
dict_table_qv = dict()

for hap in [1, 2]:
    table_qv_total = pd.read_csv(path_qv_total_result_format.format(hap = hap), names = ["seqname", "wrong", "total", "qv", "error"], sep = '\t')
    table_qv = pd.read_csv(path_qv_per_window_result_format.format(hap = hap), names = ["seqname", "wrong", "total", "qv", "error"], sep = '\t')

    table_qv.loc[np.isinf(table_qv["qv"]), "qv"] = 100
    table_qv_total.loc[np.isinf(table_qv_total["qv"]), "qv"] = 100

    table_qv["chrname"] = table_qv["seqname"].apply(lambda val: val.split(':')[-2])
    table_qv["range_start"] = table_qv["seqname"].apply(lambda val: int(val.split(':')[-1].split('-')[0]))
    table_qv["range_end"] = table_qv["seqname"].apply(lambda val: int(val.split(':')[-1].split('-')[1]))
    table_qv["range_mid"] = (table_qv["range_start"]+table_qv["range_end"]) / 2
    
    dict_table_qv_total[hap] = table_qv_total
    dict_table_qv[hap] = table_qv
# %%
list_chrname = list(map(lambda val: f"chr{val}", range(1, 23))) + ["chrX", "chrY"]
# %%
dict_chrname_to_maxlen = dict()
dict_chrname_to_qvtot = dict()

for hap in [1, 2]:
    dict_chrname_to_maxlen[hap] = dict()
    
    tmp = dict_table_qv[hap].copy()
    for chrname in list_chrname:
        if chrname not in dict_table_qv[hap]["chrname"].to_list():
            continue
        
        dict_chrname_to_maxlen[hap][chrname] = tmp[tmp["chrname"] == chrname]["range_end"].max()
    dict_chrname_to_qvtot[hap] = dict(zip(dict_table_qv_total[hap]["seqname"], dict_table_qv_total[hap]["qv"]))
# %%
plt.rcParams['font.family'] = 'Sans-serif'
plt.rcParams['font.size'] = 14

ax_height = 100
ax_height_pad = 50


total_qv_plot_ratio = 0.5
ax_width_pad_ratio = 0.125
ax_width_space_ratio = 0.2

length_divider = 1e4
xlim_padding = 2e6
qv_alert = 40

col_length_of_scatter = math.ceil(max(max(dict_chrname_to_maxlen[1].values()), max(dict_chrname_to_maxlen[2].values()))/length_divider)

fig = plt.figure(figsize=(12, 12))
row = (ax_height+ax_height_pad) * len(list_chrname) + ax_height
col = math.ceil(col_length_of_scatter*((1+total_qv_plot_ratio+ax_width_pad_ratio)*2+ax_width_space_ratio))

totalqv_width = math.floor(col_length_of_scatter*total_qv_plot_ratio)
ax_width_pad = math.floor(col_length_of_scatter*ax_width_pad_ratio)
ax_width_space = math.floor(col_length_of_scatter*ax_width_space_ratio)

# fontsize_label = 17
gsfig = GridSpec(
    row, col, 
    left=0, right=1, bottom=0,
    top=1, wspace=1, hspace=1)

dict_ax = dict()
dict_ax_bar = dict()

ax_width_start = 0
for hap in [1, 2]:
    dict_ax[hap] = dict()
    
    max_width_point = 0
    for ind, chrname in enumerate(list(filter(lambda val: val in list(dict_chrname_to_maxlen[hap].keys()), list_chrname))):
        gsrange = gsfig[
            (ax_height+ax_height_pad)*ind:(ax_height+ax_height_pad)*ind+ax_height, ax_width_start+totalqv_width+ax_width_pad:ax_width_start+totalqv_width+ax_width_pad+math.ceil(dict_chrname_to_maxlen[hap][chrname]/length_divider)
        ]
        if ax_width_start+totalqv_width+ax_width_pad+math.ceil(dict_chrname_to_maxlen[hap][chrname]/length_divider) > max_width_point:
            max_width_point = ax_width_start+totalqv_width+ax_width_pad+math.ceil(dict_chrname_to_maxlen[hap][chrname]/length_divider)
        dict_ax[hap][chrname] = fig.add_subplot(gsrange)
    
    gsrange_bar = gsfig[0:(ax_height+ax_height_pad)*(len(list(filter(lambda val: val in list(dict_chrname_to_maxlen[hap].keys()), list_chrname)))-1)+ax_height, ax_width_start:totalqv_width+ax_width_start]
    dict_ax_bar[hap] = fig.add_subplot(gsrange_bar)

    ax_width_start += (ax_width_space+max_width_point)

for hap in [1, 2]:
    list_chrname_hap = list(filter(lambda val: val in list(dict_chrname_to_maxlen[hap].keys()), list_chrname))
    
    table_qv = dict_table_qv[hap].copy()
    
    for chrname in list_chrname_hap:
        table_qv_chr = table_qv[table_qv["chrname"] == chrname]

        color = list(map(lambda val: "gray" if val > qv_alert else "r", table_qv_chr["qv"]))
        size = list(map(lambda val: 10 if val > qv_alert else 1.5, table_qv_chr["qv"]))
        dict_ax[hap][chrname].scatter(
            table_qv_chr["range_mid"], table_qv_chr["qv"],
            s = size, color = color, zorder = 5
        )
        dict_ax[hap][chrname].scatter(
            table_qv_chr[table_qv_chr["qv"]<=qv_alert]["range_mid"], table_qv_chr[table_qv_chr["qv"]<=qv_alert]["qv"],
            s = 100, color = 'r', zorder = 7, marker = '*'
        )
        dict_ax[hap][chrname].plot(
            table_qv_chr["range_mid"], table_qv_chr["qv"],
            linewidth = 0.75, color = "gray", zorder = 3
        )

        dict_ax[hap][chrname].set_yticks([0, 100])
        dict_ax[hap][chrname].set_ylim(0, 110)
        dict_ax[hap][chrname].set_xlim(0-xlim_padding, dict_chrname_to_maxlen[hap][chrname]+xlim_padding)
        dict_ax[hap][chrname].set_xticks([])
        dict_ax[hap][chrname].axhline(qv_alert, linewidth = 1, linestyle = "--", color = "gray", zorder = 1)
        # dict_ax[chrname].set_ylabel("QV")
        
        dict_ax[hap][chrname].annotate(chrname, (1+0.01 * max(dict_chrname_to_maxlen[hap].values())/dict_chrname_to_maxlen[hap][chrname], 0.5), va = "center", ha = "left", xycoords = "axes fraction", fontsize = plt.rcParams["font.size"], fontweight = "bold")
        
        dict_ax[hap][chrname].tick_params(axis = "y", labelsize = plt.rcParams["font.size"]-4)


    dict_ax_bar[hap].barh(
        list_chrname_hap[::-1], list(map(lambda name: 100 if np.isinf(dict_chrname_to_qvtot[hap][name]) else dict_chrname_to_qvtot[hap][name], list_chrname_hap))[::-1],
        color = "gray"
    )
    for ind, chrname in enumerate(list_chrname_hap[::-1]):
        dict_ax_bar[hap].annotate("inf" if dict_chrname_to_qvtot[hap][chrname] == 100 else f"{dict_chrname_to_qvtot[hap][chrname]:.2f}", (dict_chrname_to_qvtot[hap][chrname]+2, ind), va = "center", ha = "right", fontsize = plt.rcParams["font.size"]+1)
    dict_ax_bar[hap].invert_xaxis()
    dict_ax_bar[hap].set_yticks([])
    dict_ax_bar[hap].set_ylim([-0.5, len(list_chrname_hap)-0.5])
    dict_ax_bar[hap].set_xlim(150, 0)
    dict_ax_bar[hap].set_xlabel("QV per Chromosome")
    
    dict_ax_bar[hap].text(
        -0.05, 1.01, chr(ord('A')+hap-1), fontsize = plt.rcParams["font.size"]+7, fontweight = "bold", ha = "right", va = "bottom", transform = dict_ax_bar[hap].transAxes
    )
plt.show()