#%%
from pathlib import Path
import sys,gzip,math
from itertools import chain

from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import pandas as pd
import numpy as np

sys.path.append(str(Path(__file__).parents[1]))
from fasta_parser import read_fasta_as_dict

#%%
dict_hap_to_triobin = {
    1:'p',
    2:'m'
}

ax_height_gci = 75
ax_height_craq = 25
ax_height_pad = 50

length_divider = 1e4
xlim_padding = 1e6

color_qv_issue = "darkorange"
color_gci_depth = "cornflowerblue"
color_gci_issue = "crimson"
color_craq_cre_issue = "mediumblue"
color_craq_cse_issue = "magenta"

depth_window = 50000

path_gci_depth_format = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/GCI_Evaluation/GCI_Result/Hap{hap}/GCI_Result.hap{hap}.HiFi_ONT.Winnowmap2_Minimap2.triobin_{triobin}0a_two_type.depth.mean_50000.txt"
path_gci_issue_format = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/GCI_Evaluation/GCI_Result/Hap{hap}/GCI_Result.hap{hap}.HiFi_ONT.Winnowmap2_Minimap2.triobin_{triobin}0a_two_type.0.depth.bed"
path_craq_cre_format = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/CRAQ_Evaluation/CRAQ_Result/LongRead_Merged_All.ShortRead/Hap{hap}/runAQI_out/locER_out/out_final.CRE.bed"
path_craq_cse_format = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/CRAQ_Evaluation/CRAQ_Result/LongRead_Merged_All.ShortRead/Hap{hap}/runAQI_out/strER_out/out_final.CSE.bed"
path_qv_issue_format = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Merqury_Evaluation/Merqury_QV_per_window/PolishedASM.KOREF_S1.Trio_Specific.UpperCase.Polish_R8.SV_SNV.hap{hap}.QV_per_window.k21.chunked.10000.qv"
path_fasta_format = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Final_Assembly/FinalASM.KOREF_S1.Trio_Specific.hap{hap}.Reordered.fasta"


def read_qv_issue(path_qv_issue):
    dict_chr_to_qv = dict()
    
    with open(path_qv_issue, 'r') as fr:
        for line in fr:
            values = line.strip('\n').split()
            regionname = values[0]
            qv_local = values[3]
            qv_local = None if qv_local == "inf" else float(qv_local)

            chrname = regionname.split(':')[-2]
            start, end = regionname.split(':')[-1].split('-')
            
            if dict_chr_to_qv.get(chrname) == None:
                dict_chr_to_qv[chrname] = dict()
            if qv_local != None:
                dict_chr_to_qv[chrname][(int(start)+1, int(end))] = qv_local
    return dict_chr_to_qv

def draw_issue_region_module(path_fasta, path_gci_depth, path_gci_issue, path_craq_cre, path_craq_cse, path_qv_issue, ax, draw_legend = True):
    dict_fasta = read_fasta_as_dict(path_fasta)
    list_chromosomes = list(dict_fasta.keys())

    dict_gci_windowed_depth = dict()
    with open(path_gci_depth, "r") as fr:
        chrname = None
        
        for line in fr:
            line = line.strip('\n')
            if line.startswith('>'):
                chrname = line[1:]
                print(chrname)
                assert dict_gci_windowed_depth.get(chrname) == None
                dict_gci_windowed_depth[chrname] = list()
            else:
                val = float(line)
                dict_gci_windowed_depth[chrname].append(val)

    table_gci_issue = pd.read_csv(path_gci_issue, sep = '\t', names = ["chr", "start", "end"])
    table_craq_cre_issue = pd.read_csv(path_craq_cre, sep = '\t', names = ["chr", "start", "end", "name", "type"])
    table_craq_cse_issue = pd.read_csv(path_craq_cse, sep = '\t', names = ["chr", "start", "end", "name", "type"])
    
    dict_qv_issue = read_qv_issue(path_qv_issue)

    dict_chrname_to_maxlen = dict()
    for chrname, sequence in dict_fasta.items():
        dict_chrname_to_maxlen[chrname] = len(sequence)

    
    table_craq_cre_issue["Region"] = (table_craq_cre_issue["end"]-table_craq_cre_issue["start"]).apply(lambda val: max(1, val))
    table_craq_cse_issue["Region"] = (table_craq_cse_issue["end"]-table_craq_cse_issue["start"]).apply(lambda val: max(1, val))

    col_length_of_scatter = math.ceil(max(dict_chrname_to_maxlen.values())/length_divider)

    # fig = plt.figure(figsize=(12, 15))
    ax.axis("off")
    row = (ax_height_gci+ax_height_craq+ax_height_pad) * (len(list_chromosomes)-1) + ax_height_gci+ax_height_craq
    col = math.ceil(col_length_of_scatter)

    gsfig = GridSpec(
        row, col, 
        left=0, right=1, bottom=0,
        top=1, wspace=1, hspace=1)

    dict_ax_gci = dict()
    dict_ax_craq = dict()

    gci_ylim_max = np.quantile(list(chain(*dict_gci_windowed_depth.values())), 0.999)

    for ind, chrname in enumerate(list_chromosomes):
        gsrange_gci = gsfig[(
            ax_height_gci+ax_height_craq+ax_height_pad)*ind:(ax_height_gci+ax_height_craq+ax_height_pad)*ind+ax_height_gci, 
            0:0+math.ceil((dict_chrname_to_maxlen[chrname]+2*xlim_padding)/length_divider
        )]
        dict_ax_gci[chrname] = ax.inset_axes([gsrange_gci.colspan.start/col, 1-gsrange_gci.rowspan.stop/row, (gsrange_gci.colspan.stop-gsrange_gci.colspan.start)/col, (gsrange_gci.rowspan.stop-gsrange_gci.rowspan.start)/row])
        
        gsrange_craq = gsfig[(
            ax_height_gci+ax_height_craq+ax_height_pad)*ind+ax_height_gci:(ax_height_gci+ax_height_craq+ax_height_pad)*ind+ax_height_gci+ax_height_craq, 
            0:0+math.ceil((dict_chrname_to_maxlen[chrname]+2*xlim_padding)/length_divider
        )]
        dict_ax_craq[chrname] = ax.inset_axes([gsrange_craq.colspan.start/col, 1-gsrange_craq.rowspan.stop/row, (gsrange_craq.colspan.stop-gsrange_craq.colspan.start)/col, (gsrange_craq.rowspan.stop-gsrange_craq.rowspan.start)/row])

    ########## Plot GCI Depth
    for ind, chrname in enumerate(list_chromosomes):
        max_ntimes_window = math.floor(dict_chrname_to_maxlen[chrname] / depth_window)
        
        list_average_window_pos = [0] + list(depth_window/2 + np.arange(0, max_ntimes_window) * depth_window)
        list_average_window_pos.append((max_ntimes_window * depth_window + dict_chrname_to_maxlen[chrname])/2)
        list_average_window_pos.append(dict_chrname_to_maxlen[chrname])
        
        list_average_depth = [0] + dict_gci_windowed_depth[chrname] + [0]
        
        combined_depth = dict_ax_gci[chrname].fill_between(
            list_average_window_pos, list_average_depth, color = color_gci_depth, alpha = 0.5, zorder = 1
        )
        
        for _, row in table_gci_issue[table_gci_issue["chr"] == chrname].iterrows():
            range_start = row["start"]
            range_end = row["end"]
            
            if range_start == dict_chrname_to_maxlen[chrname] or range_end == 1:
                continue
            
            if range_start == range_end:
                dict_ax_gci[chrname].axvline(range_start, linewidth = 1, color = color_gci_issue, alpha = 0.5, zorder = 3)
            else:
                gci_issue = dict_ax_gci[chrname].axvspan(
                    range_start, range_end, color = color_gci_issue, alpha = 0.5, zorder = 3
                )
        
        dict_ax_gci[chrname].set_xlim(-xlim_padding, dict_chrname_to_maxlen[chrname]+xlim_padding)
        dict_ax_gci[chrname].set_ylim(0, gci_ylim_max)
        dict_ax_gci[chrname].set_xticks([])
        
        dict_ax_gci[chrname].annotate(chrname, (1+0.01 * max(dict_chrname_to_maxlen.values())/dict_chrname_to_maxlen[chrname], 0.5*ax_height_craq/(ax_height_gci-ax_height_craq)), va = "center", ha = "left", xycoords = "axes fraction", fontsize = plt.rcParams["font.size"]+5, fontweight = "bold")

    ########## Plot CRAQ
    for ind, chrname in enumerate(list_chromosomes):
        
        for _, row in table_craq_cre_issue[table_craq_cre_issue["chr"] == chrname].iterrows():
            range_start = row["start"]
            range_end = row["end"]
            
            if range_start == dict_chrname_to_maxlen[chrname] or range_end == 1:
                continue
            
            if range_start == range_end:
                dict_ax_craq[chrname].axvline(range_start, linewidth = 1, color = color_craq_cre_issue, alpha = 1, zorder = 3)
            else:
                craq_cre = dict_ax_craq[chrname].axvspan(
                    range_start, range_end, color = color_craq_cre_issue, alpha = 0.5, zorder = 3
                )
        for _, row in table_craq_cse_issue[table_craq_cse_issue["chr"] == chrname].iterrows():
            range_start = row["start"]
            range_end = row["end"]
            
            if range_start == dict_chrname_to_maxlen[chrname] or range_end == 1:
                continue
            
            if range_start == range_end:
                dict_ax_craq[chrname].axvline(range_start, linewidth = 1, color = color_craq_cse_issue, alpha = 1, zorder = 3)
            else:
                craq_cse = dict_ax_craq[chrname].axvspan(
                    range_start, range_end, color = color_craq_cse_issue, alpha = 0.5, zorder = 3
                )
        
        dict_ax_craq[chrname].set_xlim(-xlim_padding, dict_chrname_to_maxlen[chrname]+xlim_padding)
        dict_ax_craq[chrname].set_xticks([])
        dict_ax_craq[chrname].set_yticks([])
    
    ########## Plot QV Not Inf
    for ind, chrname in enumerate(list_chromosomes):
        for (start, end), qv in dict_qv_issue[chrname].items():
            qvissue = dict_ax_craq[chrname].axvspan(
                start, end, color = color_qv_issue, alpha = 0.5, zorder = 3
            )

    if draw_legend:
        legend_gci = ax.legend(handles = [combined_depth, gci_issue], labels = ["Combined Depth (ONT+HiFi)", "Low-Depth Issue (<1)"], loc = "upper center", bbox_to_anchor = (0.85, 0.30), bbox_transform = ax.transAxes, title = "GCI Report", fontsize = plt.rcParams["font.size"]+2, title_fontproperties = {"weight":"bold", "size":plt.rcParams["font.size"]+4})

        legend_craq = ax.legend(handles = [craq_cse, craq_cre, qvissue], labels = ["Large Structural Breakage (CSE)", "Regional Error (CRE)", "QV Issue (Not Infinite)"], loc = "upper center", bbox_to_anchor = (0.85, 0.20), bbox_transform = ax.transAxes, title = "CRAQ & Merqury Report", fontsize = plt.rcParams["font.size"]+2, title_fontproperties = {"weight":"bold", "size":plt.rcParams["font.size"]+4})

        # legend_qv = ax.legend(handles = [qvissue], labels = ["QV Issue (Not Infinite)"], loc = "upper center", bbox_to_anchor = (0.85, 0.15), bbox_transform = ax.transAxes, title = "Merqury Report", fontsize = plt.rcParams["font.size"]+2, title_fontproperties = {"weight":"bold", "size":plt.rcParams["font.size"]+4})

        ax.add_artist(legend_gci)
        # ax.add_artist(legend_craq)
# %%
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams["font.size"] = 11

fig = plt.figure(figsize=(15, 15))
row = 100
col = 100
fontsize_label = 30
gsfig = GridSpec(
    row, col, 
    left=0, right=1, bottom=0,
    top=1, wspace=1, hspace=1)


gs1 = gsfig[0:100, 0:45] 
ax1 = fig.add_subplot(gs1)

gs2 = gsfig[0:100, 55:100]
ax2 = fig.add_subplot(gs2)

hap_ax1 = 1
draw_issue_region_module(
    path_fasta_format.format(hap = hap_ax1), path_gci_depth_format.format(hap = hap_ax1, triobin =dict_hap_to_triobin[hap_ax1]), path_gci_issue_format.format(hap = hap_ax1, triobin =dict_hap_to_triobin[hap_ax1]), path_craq_cre_format.format(hap = hap_ax1), path_craq_cse_format.format(hap = hap_ax1), path_qv_issue_format.format(hap = hap_ax1), ax1,
    draw_legend = False
)
ax1.text(-0.03, 1.001, 'A', fontdict={"size":fontsize_label, "weight":"bold"}, ha = "right", va = "bottom")

hap_ax2 = 2
draw_issue_region_module(
    path_fasta_format.format(hap = hap_ax2), path_gci_depth_format.format(hap = hap_ax2, triobin =dict_hap_to_triobin[hap_ax2]), path_gci_issue_format.format(hap = hap_ax2, triobin =dict_hap_to_triobin[hap_ax2]), path_craq_cre_format.format(hap = hap_ax2), path_craq_cse_format.format(hap = hap_ax2), path_qv_issue_format.format(hap = hap_ax2), ax2
)
ax2.text(-0.03, 1.001, 'B', fontdict={"size":fontsize_label, "weight":"bold"}, ha = "right", va = "bottom")

plt.show()
# %%
