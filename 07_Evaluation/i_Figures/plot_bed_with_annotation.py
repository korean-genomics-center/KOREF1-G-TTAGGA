#%%
import math
from collections import Counter

import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec

# %%
path_bed = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Repeat_Annotation/Segmental_Duplication/Biser/Hap1/FinalASM.KOREF_S1.Trio_Specific.hap1.Reordered.masked.SD.check_overlap.filt_overlap_50.err_lt30.merged.bed"

path_fai = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Final_Assembly/FinalASM.KOREF_S1.Trio_Specific.hap1.Reordered.fasta.fai"

#%%
draw_annot = False

col_annot = 4
name_annot = "Repeat Type"

annot_name_conv = {
    "DNA" : "DNA Elements",
    "rRNA" : "Small RNA",
    "snRNA" : "Small RNA",
    "srpRNA" : "Small RNA",
    "scRNA" : "Small RNA",
    "tRNA" : "Small RNA"
}
#%%
dict_chr_to_len = dict()

with open(path_fai, 'r') as fr:
    for line in fr:
        values = line.strip('\n').split()
        dict_chr_to_len[values[0]] = int(values[1])

with open(path_bed, 'r') as fr:
    while 1:
        fstline = fr.readline()
        if '\t' in fstline: break
    n_elements = len(fstline.split('\t'))

table_bed = pd.read_csv(path_bed, sep = '\t', names = ["chr", "start", "end"] + list(range(1, n_elements+1-3)))
# %%
if draw_annot:
    table_bed["annot"] = table_bed[col_annot].apply(lambda val: annot_name_conv.get(val, val))

    color_palette = sns.color_palette("Paired", 12)
    dict_palette = dict(zip(table_bed["annot"].unique(), color_palette))
# %%
plt.rcParams['font.family'] = 'Sans-serif'
plt.rcParams['font.size'] = 14

plt.figure(figsize = (4, 3))
plt.hist(table_bed["end"] - table_bed["start"], bins = 100, color = "gray")
plt.xlabel("Length (bp)")
plt.ylabel("Count")
plt.yscale("log")
plt.show()

#%%
plt.rcParams['font.family'] = 'Sans-serif'
plt.rcParams['font.size'] = 14

ax_height = 50
ax_height_pad = 50

length_divider = 1e4
xlim_padding = 1e5

min_length_draw = 100

col_length_of_scatter = math.ceil(max(dict_chr_to_len.values())/length_divider)


fig = plt.figure(figsize=(8, 7))

row = (ax_height+ax_height_pad) * len(dict_chr_to_len) + ax_height
col = math.ceil(col_length_of_scatter)

gsfig = GridSpec(
    row, col, 
    left=0, right=1, bottom=0,
    top=1, wspace=1, hspace=1)

dict_ax = dict()

ax_width_start = 0

    
for ind, chrname in enumerate(list(dict_chr_to_len.keys())):
    gsrange = gsfig[
        (ax_height+ax_height_pad)*ind:(ax_height+ax_height_pad)*ind+ax_height, ax_width_start:ax_width_start+math.ceil((dict_chr_to_len[chrname]+xlim_padding*2)/length_divider)
    ]
    dict_ax[chrname] = fig.add_subplot(gsrange)


for ind, chrname in enumerate(list(dict_chr_to_len.keys())):
    table_bed_chr = table_bed[table_bed["chr"] == chrname]

    if draw_annot:
        for annot in dict_palette.keys():
            table_bed_chr_annot = table_bed_chr[table_bed_chr["annot"] == annot]
            for _, row in table_bed_chr_annot.iterrows():
                repeatstart, repeatend = sorted([row["start"], row["end"]])
                # if (repeatend-repeatstart) < min_length_draw:
                #     continue

                if ind == len(dict_chr_to_len)-1:
                    dict_ax[chrname].barh([0], [repeatend-repeatstart], left = [repeatstart], color = dict_palette[annot], label = annot)
                else:
                    dict_ax[chrname].barh([0], [repeatend-repeatstart], left = [repeatstart], color = dict_palette[annot])
                # break
    else:
        for _, row in table_bed_chr.iterrows():
            repeatstart, repeatend = sorted([row["start"], row["end"]])

            dict_ax[chrname].barh([0], [repeatend-repeatstart], left = [repeatstart], color = "gray")

    dict_ax[chrname].set_xlim(0-xlim_padding, dict_chr_to_len[chrname]+xlim_padding)
    dict_ax[chrname].set_yticks([])
    dict_ax[chrname].set_xticks([])
    dict_ax[chrname].annotate(chrname, (-0.01 * max(dict_chr_to_len.values())/dict_chr_to_len[chrname], 0.5), va = "center", ha = "right", xycoords = "axes fraction", fontsize = plt.rcParams["font.size"], fontweight = "bold")

if draw_annot:
    plt.legend(loc = "lower right", bbox_to_anchor = (0.98, 0.02), bbox_transform = fig.transFigure)
plt.show()
# %%
