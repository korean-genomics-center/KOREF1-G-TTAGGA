#%%
from pathlib import Path
import sys,gzip,math
from itertools import chain

from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import pandas as pd
import numpy as np

sys.path.append(str(Path(__file__).parents[2] / "common"))
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

#%%
dict_qv = {
    "KOREF1-G-TTAGGA (Pat)": 81.19,
    "KOREF1-G-TTAGGA (Mat)": 79.03,
    "CHM13": 73.94,
    "CN1 (Pat)": 59.36,
    "CN1 (Mat)": 60.10,
    "T2T-YAO (Pat)": 72.28,
    "T2T-YAO (Mat)": 70.49,
    "I002C (Pat)": 83.08,
    "I002C (Mat)": 82.05
}

dict_chr_order = dict(zip(
    list(map(lambda val: f"chr{val}", range(1, 23))) + ["chrX", "chrY"],
    list(range(24))
))

path_qv_per_chr_pat = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Merqury_Evaluation/Merqury_QV/PolishedASM.KOREF_S1.Trio_Specific.UpperCase.Polish_R8.SV_SNV.hap1.k21.perseq.qv"
path_qv_per_chr_mat = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Merqury_Evaluation/Merqury_QV/PolishedASM.KOREF_S1.Trio_Specific.UpperCase.Polish_R8.SV_SNV.hap2.k21.perseq.qv"

dict_path_fasta = {
    "KOREF1-G-TTAGGA (Pat)": "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Final_Assembly/FinalASM.KOREF_S1.Trio_Specific.hap1.Reordered.fasta",
    "KOREF1-G-TTAGGA (Mat)": "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Final_Assembly/FinalASM.KOREF_S1.Trio_Specific.hap2.Reordered.fasta",
    "CHM13": "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/CHM13/chm13v2.0.fa",
    "CN1 (Pat)": "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/CN1/CN1_pat.v1.0.1.fasta",
    "CN1 (Mat)": "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/CN1/CN1_mat.v1.0.1.fasta",
    "T2T-YAO (Pat)": "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/T2T_YAO/T2T_Yao_pat_v1.1.fasta",
    "T2T-YAO (Mat)": "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/T2T_YAO/T2T_Yao_mat_v1.1.fasta",
    "I002C (Pat)": "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/I002C/I002C_Paternal_v0.7.fasta",
    "I002C (Mat)": "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/I002C/I002C_Maternal_v0.7.fasta"
}

#%%
dict_asmname_to_fasta = dict()

for asmname, path in dict_path_fasta.items():
    dict_asmname_to_fasta[asmname] = read_fasta_as_dict(path)
# %%
def get_total_length(dict_fasta):
    return sum(list(map(len, dict_fasta.values())))

def get_ng50(dict_fasta):
    lengths_chr = sorted(list(map(len, dict_fasta.values())))
    totlen = sum(lengths_chr)
    
    cumsum = 0
    ng50 = None
    for length in lengths_chr:
        cumsum += length
        if cumsum >= totlen/2:
            ng50 = length
            break
    return ng50

#%%
dict_asm_to_totlen = dict()
dict_asm_to_ng50 = dict()

for asmname, dict_fasta in dict_asmname_to_fasta.items():
    dict_asm_to_totlen[asmname] = get_total_length(dict_fasta)
    dict_asm_to_ng50[asmname] = get_ng50(dict_fasta)
# %%
table_qv_per_seq_pat = pd.read_csv(path_qv_per_chr_pat, delim_whitespace = True, names = ["chr", "nerr", "tot", "qv", "errrate"])
table_qv_per_seq_mat = pd.read_csv(path_qv_per_chr_mat, delim_whitespace = True, names = ["chr", "nerr", "tot", "qv", "errrate"])

table_qv_per_seq_pat.loc[np.isposinf(table_qv_per_seq_pat["qv"]), "qv"] = 100
table_qv_per_seq_mat.loc[np.isposinf(table_qv_per_seq_mat["qv"]), "qv"] = 100

table_qv_per_seq_pat["chrorder"] = table_qv_per_seq_pat["chr"].apply(dict_chr_order.__getitem__)
table_qv_per_seq_mat["chrorder"] = table_qv_per_seq_mat["chr"].apply(dict_chr_order.__getitem__)

table_qv_per_seq_pat = table_qv_per_seq_pat.sort_values(by = ["chrorder"]).reset_index(drop = True)
table_qv_per_seq_mat = table_qv_per_seq_mat.sort_values(by = ["chrorder"]).reset_index(drop = True)
#%%
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
plt.rcParams["font.size"] = 14
fontsize_label = 20

fig = plt.figure(figsize=(12, 6))
row = 100
col = 101
gsfig = GridSpec(
    row, col, 
    left=0, right=1, bottom=0,
    top=1, wspace=1, hspace=1)


gs_qv1 = gsfig[0:40, 0:47] 
ax_qv1 = fig.add_subplot(gs_qv1)

gs_qv2 = gsfig[0:40, 101-47:101]
ax_qv2 = fig.add_subplot(gs_qv2)

hspace_layer2 = 7
width_layer2 = 29

gs_totlen = gsfig[57:100, 0:width_layer2]
ax_totlen = fig.add_subplot(gs_totlen)

gs_ng50 = gsfig[57:100, width_layer2+hspace_layer2:width_layer2*2+hspace_layer2]
ax_ng50 = fig.add_subplot(gs_ng50)

gs_qvcomp = gsfig[57:100, width_layer2*2+hspace_layer2*2:width_layer2*3+hspace_layer2*2]
ax_qvcomp = fig.add_subplot(gs_qvcomp)


label_x = -0.03
label_y = 0.01

########### QV per chr (Pat)
ax_qv1.bar(
    table_qv_per_seq_pat["chr"],
    table_qv_per_seq_pat["qv"],
    color = "gray",
    zorder = 3
)
for ind, row in table_qv_per_seq_pat.iterrows():
    ax_qv1.annotate(f"{row['qv']:.1f}" if row["qv"] < 100 else "", (ind, row["qv"]+2), rotation = 90, ha = "center", va = "bottom")

ax_qv1.set_xticks(range(table_qv_per_seq_pat.shape[0]), table_qv_per_seq_pat["chr"], rotation = -90, rotation_mode = "anchor", ha = "left", va = "center")
ax_qv1.set_xlim(-0.5, table_qv_per_seq_pat.shape[0]-0.5)
ax_qv1.set_ylim(0, 120)
ax_qv1.set_yticks(np.linspace(0, 100, 6))
ax_qv1.grid(axis = 'y', zorder = 1, color = "lightgray", linewidth = 1, linestyle = "--")
ax_qv1.set_ylabel("QV")
ax_qv1.set_title("Paternal Assembly", fontweight = "bold")
ax_qv1.text(
    label_x, 1+label_y, 'A', fontsize = fontsize_label, fontweight = "bold", transform = ax_qv1.transAxes, ha = "right", va = "bottom"
)

########### QV per chr (Mat)
ax_qv2.bar(
    table_qv_per_seq_mat["chr"],
    table_qv_per_seq_mat["qv"],
    color = "gray",
    zorder = 3
)
for ind, row in table_qv_per_seq_mat.iterrows():
    ax_qv2.annotate(f"{row['qv']:.1f}" if row["qv"] < 100 else "", (ind, row["qv"]+2), rotation = 90, ha = "center", va = "bottom")

ax_qv2.set_xticks(range(table_qv_per_seq_mat.shape[0]), table_qv_per_seq_mat["chr"], rotation = -90, rotation_mode = "anchor", ha = "left", va = "center")
ax_qv2.set_xlim(-0.5, table_qv_per_seq_mat.shape[0]-0.5)
ax_qv2.set_ylim(0, 120)
ax_qv2.set_yticks(np.linspace(0, 100, 6))
ax_qv2.grid(axis = 'y', zorder = 1, color = "lightgray", linewidth = 1, linestyle = "--")
ax_qv2.set_ylabel("QV")
ax_qv2.set_title("Maternal Assembly", fontweight = "bold")
ax_qv2.text(
    label_x * ((gs_qv1.colspan.stop-gs_qv1.colspan.start) / (gs_qv2.colspan.stop-gs_qv2.colspan.start)), 1+label_y*((gs_qv1.rowspan.stop-gs_qv1.rowspan.start) / (gs_qv2.rowspan.stop-gs_qv2.rowspan.start)), 'B', fontsize = fontsize_label, fontweight = "bold", transform = ax_qv2.transAxes, ha = "right", va = "bottom"
)

########### TotLength
ax_totlen.bar(
    dict_asm_to_totlen.keys(),
    dict_asm_to_totlen.values(),
    color = "gray",
    zorder = 3
)
for ind, asmlen in enumerate(dict_asm_to_totlen.values()):
    ax_totlen.annotate(
        f"{asmlen/1e9:.2f}G",
        (ind, asmlen+0.1e9), rotation = 90, ha = "center", va = "bottom"
    )
    
ax_totlen.set_yticks(np.linspace(0, 3e9, 4))
ax_totlen.grid(axis = 'y', zorder = 1, color = "lightgray", linewidth = 1, linestyle = "--")
ax_totlen.set_xticks(range(len(dict_asm_to_totlen)), dict_asm_to_totlen.keys(), rotation = -45, rotation_mode = "anchor", ha = "left", va = "center")
ax_totlen.set_ylim(0, 4.2e9)
ax_totlen.set_xlim(-0.5, len(dict_asm_to_totlen)-0.5)
ax_totlen.set_ylabel("size (bp)")
ax_totlen.set_title("Assembly Size", fontdict = {"weight":"bold"})
ax_totlen.text(
    label_x * ((gs_qv1.colspan.stop-gs_qv1.colspan.start) / (gs_totlen.colspan.stop-gs_totlen.colspan.start)), 1+label_y*((gs_qv1.rowspan.stop-gs_qv1.rowspan.start) / (gs_totlen.rowspan.stop-gs_totlen.rowspan.start)), 'C', fontsize = fontsize_label, fontweight = "bold", transform = ax_totlen.transAxes, ha = "right", va = "bottom"
)

########### NG50
ax_ng50.bar(
    dict_asm_to_ng50.keys(),
    dict_asm_to_ng50.values(),
    color = "gray",
    zorder = 3
)
for ind, ng50 in enumerate(dict_asm_to_ng50.values()):
    ax_ng50.annotate(
        f"{ng50/1e6:.1f}M",
        (ind, ng50+0.05e8), rotation = 90, ha = "center", va = "bottom"
    )
    
ax_ng50.set_yticks(np.linspace(0, 1.5e8, 4))
ax_ng50.grid(axis = 'y', zorder = 1, color = "lightgray", linewidth = 1, linestyle = "--")
ax_ng50.set_xticks(range(len(dict_asm_to_totlen)), dict_asm_to_totlen.keys(), rotation = -45, rotation_mode = "anchor", ha = "left", va = "center")
ax_ng50.set_ylim(0, 2.25e8)
ax_ng50.set_xlim(-0.5, len(dict_asm_to_totlen)-0.5)
ax_ng50.set_ylabel("size (bp)")
ax_ng50.set_title("NG50", fontdict = {"weight":"bold"})
ax_ng50.text(
    label_x * ((gs_qv1.colspan.stop-gs_qv1.colspan.start) / (gs_ng50.colspan.stop-gs_ng50.colspan.start)), 1+label_y*((gs_qv1.rowspan.stop-gs_qv1.rowspan.start) / (gs_ng50.rowspan.stop-gs_ng50.rowspan.start)), 'D', fontsize = fontsize_label, fontweight = "bold", transform = ax_ng50.transAxes, ha = "right", va = "bottom"
)

########### QV
ax_qvcomp.bar(
    dict_qv.keys(),
    dict_qv.values(),
    color = "gray",
    zorder = 3
)
for ind, qv in enumerate(dict_qv.values()):
    ax_qvcomp.annotate(
        f"{qv:.1f}",
        (ind, qv+2), rotation = 90, ha = "center", va = "bottom"
    )
    
ax_qvcomp.set_yticks(np.linspace(0, 80, 5))
ax_qvcomp.grid(axis = 'y', zorder = 1, color = "lightgray", linewidth = 1, linestyle = "--")
ax_qvcomp.set_xticks(range(len(dict_asm_to_totlen)), dict_asm_to_totlen.keys(), rotation = -45, rotation_mode = "anchor", ha = "left", va = "center")
ax_qvcomp.set_ylim(0, 105)
ax_qvcomp.set_xlim(-0.5, len(dict_asm_to_totlen)-0.5)
ax_qvcomp.set_ylabel("QV")
ax_qvcomp.set_title("QV", fontdict = {"weight":"bold"})
ax_qvcomp.text(
    label_x * ((gs_qv1.colspan.stop-gs_qv1.colspan.start) / (gs_qvcomp.colspan.stop-gs_qvcomp.colspan.start)), 1+label_y*((gs_qv1.rowspan.stop-gs_qv1.rowspan.start) / (gs_qvcomp.rowspan.stop-gs_qvcomp.rowspan.start)), 'E', fontsize = fontsize_label, fontweight = "bold", transform = ax_qvcomp.transAxes, ha = "right", va = "bottom"
)
plt.show()
# %%
