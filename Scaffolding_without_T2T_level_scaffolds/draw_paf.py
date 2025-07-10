#%%
import pandas as pd
import os
import numpy as np
from collections import Counter

path_paf = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/Round6.TGS_Gapcloser_Modified_V2.Unplaced_Scf.HiFi_L1000_QV20_parental.Hifiasm_BinnedASM_pctg.Flye_HiFi_BinnedASM.Flye_ONT_BinnedASM_asm5P_hard.Manual_Scaffolding/Prealignment_for_Candidate_Filtering/Hap1/Align/Align.REF_KOREFDraftASM_Hap1_Split.QRY_KOREFDraftASM_Hap1_Split.asm5P.paf"

paf_filt_len = 300
paf_filt_quality = 0
paf_filt_identity = 0.7

#%%
table_paf = pd.read_csv(path_paf, sep = '\t', names = list(range(24)))
# %%
col_qry = 0
col_qry_len = 1
col_qry_start = 2
col_qry_end = 3
col_strand = 4
col_ref = 5
col_ref_len = 6
col_ref_start = 7
col_ref_end = 8
col_match_len = 9
col_align_len = 10
col_quality = 11
#%%
table_paf["identity"] = table_paf[col_match_len] / table_paf[col_align_len]

# %%
dict_qry_to_ref_match = dict()

for qryname in table_paf[col_qry].unique():
    table_paf_qry = table_paf[table_paf[col_qry] == qryname]
    dict_qry_to_ref_match[qryname] = dict()
    for refname in table_paf_qry[col_ref].unique():
        table_paf_qry_to_ref = table_paf_qry[table_paf_qry[col_ref] == refname]

        reflen = table_paf_qry_to_ref.loc[:, col_ref_len].to_list()[0]
        all_matched_len = table_paf_qry_to_ref.loc[:, col_match_len].sum()
        all_align_len = table_paf_qry_to_ref.loc[:, col_align_len].sum()

        strand_plus_len = table_paf_qry_to_ref[table_paf_qry_to_ref[col_strand] == '+'].loc[:, col_match_len].sum()
        strand_minus_len = table_paf_qry_to_ref[table_paf_qry_to_ref[col_strand] == '-'].loc[:, col_match_len].sum()

        if strand_plus_len > strand_minus_len:
            strand = '+'
            strand_ratio = strand_plus_len / (strand_plus_len + strand_minus_len)
        else:
            strand = '-'
            strand_ratio = strand_minus_len / (strand_plus_len + strand_minus_len)
        dict_qry_to_ref_match[qryname][refname] = {
            "matched_ratio" : all_matched_len/all_align_len,
            "strand" : strand,
            "strand_ratio" : strand_ratio
        }
# %%
table_paf_draw = table_paf[table_paf[col_align_len] >= paf_filt_len]
table_paf_draw = table_paf_draw[table_paf_draw[col_quality] >= paf_filt_quality]
table_paf_draw = table_paf_draw[table_paf_draw["identity"] >= paf_filt_identity]

# table_paf_draw = table_paf_draw[table_paf_draw[7] < 5000000]
# %%
table_paf_draw_qry_dropdup = table_paf_draw.drop_duplicates(subset = [col_qry])
dict_scf_to_scf_length = dict(zip(
    table_paf_draw_qry_dropdup[col_qry],
    table_paf_draw_qry_dropdup[col_qry_len]
))
table_paf_draw_ref_dropdup = table_paf_draw.drop_duplicates(subset = [col_ref])
dict_ref_to_ref_length = dict(zip(
    table_paf_draw_ref_dropdup[col_ref],
    table_paf_draw_ref_dropdup[col_ref_len]
))
# %%
list_scf_order = sorted(dict_scf_to_scf_length.keys(), key = dict_scf_to_scf_length.__getitem__)
list_ref_order = sorted(dict_ref_to_ref_length.keys(), key = dict_ref_to_ref_length.__getitem__)[::-1]
#%%
from matplotlib import pyplot as plt

strand_plus_color = "salmon"
strand_minus_color = "cornflowerblue"

# list_ref_order = ["chr11", "chr15"]
# list_ref_order = list(map(lambda val: f"chr{val}", range(1, 23))) + ["chrX", "chrY"]
# list_scf_order = (list(map(lambda val: f"chr{val}", range(1, 23))) + ["chrX", "chrY"])[::-1]

# list_ref_order = ["chrX", "chrY"]
# list_ref_order = ["chr13", "chr14", "chr15", "chr21", "chr22"]

equal_scf_ratio = True
equal_ref_ratio = False

plt.rcParams["font.size"] = 20
fig, axes = plt.subplots(
    len(list_scf_order),len(list_ref_order),
    figsize=(15,15), 
    # figsize=(10,10), 
    gridspec_kw={
        'height_ratios': [1] * len(list_scf_order) if equal_scf_ratio else list(map(dict_scf_to_scf_length.__getitem__, list_scf_order)),
        "width_ratios":[1] * len(list_ref_order) if equal_ref_ratio else list(map(dict_ref_to_ref_length.__getitem__, list_ref_order))}
    # gridspec_kw={
    #     'height_ratios': list(map(dict_ref_to_ref_length.__getitem__, list_ref_order))[::-1],
    #     "width_ratios":list(map(dict_ref_to_ref_length.__getitem__, list_ref_order))}
)
plt.subplots_adjust(wspace = 0.1, hspace = 0)

for ind_scf, scf in enumerate(list_scf_order):
    for ind_ref, ref in enumerate(list_ref_order):
        if ind_ref == 0:
            axes[ind_scf, ind_ref].set_ylabel(scf, rotation = 0, ha = "right", va = 'center')
        if ind_scf == len(list_scf_order)-1:
            axes[ind_scf, ind_ref].set_xlabel(ref, rotation = 90)
        axes[ind_scf, ind_ref].set_xticks([])
        axes[ind_scf, ind_ref].set_yticks([])

        table_paf_draw_of_this_comb = table_paf_draw[
            np.logical_and(
                table_paf_draw[col_ref] == ref,
                table_paf_draw[col_qry] == scf
            )
        ]
        for _, row in table_paf_draw_of_this_comb.iterrows():
            scf_start = row[col_qry_start]
            scf_end = row[col_qry_end]
            scf_size = row[col_qry_len]
            ref_start = row[col_ref_start]
            ref_end = row[col_ref_end]
            ref_size = row[col_ref_len]
            strand = row[col_strand]
            assert strand in {'+', '-'}, row

            axes[ind_scf, ind_ref].plot(
                [ref_start, ref_end],
                [scf_start, scf_end] if strand == '+' else [scf_end, scf_start],
                linewidth = 4,
                color = strand_plus_color if strand == '+' else strand_minus_color,
                clip_on = False
            )
            axes[ind_scf, ind_ref].scatter(
                [ref_start, ref_end],
                [scf_start, scf_end] if strand == '+' else [scf_end, scf_start],
                s = 15,
                color = strand_plus_color if strand == '+' else strand_minus_color,
                clip_on = False
            )
            axes[ind_scf, ind_ref].set_xlim(0, ref_size)
            axes[ind_scf, ind_ref].set_ylim(0, scf_size)
plt.show()
# %%
from matplotlib import pyplot as plt

strand_plus_color = "salmon"
strand_minus_color = "cornflowerblue"
###### Draw Single
ref = 3
scf = 4

plt.rcParams["font.size"] = 14
fig, ax = plt.subplots(1, 1, figsize = (4, 4))
table_paf_draw_of_this_comb = table_paf_draw[
    np.logical_and(
        table_paf_draw[col_ref] == ref,
        table_paf_draw[col_qry] == scf
    )
]
for _, row in table_paf_draw_of_this_comb.iterrows():
    scf_start = row[col_qry_start]
    scf_end = row[col_qry_end]
    scf_size = row[col_qry_len]
    ref_start = row[col_ref_start]
    ref_end = row[col_ref_end]
    ref_size = row[col_ref_len]
    strand = row[col_strand]
    assert strand in {'+', '-'}, row

    ax.plot(
        [ref_start, ref_end],
        [scf_start, scf_end] if strand == '+' else [scf_end, scf_start],
        linewidth = 3,
        color = strand_plus_color if strand == '+' else strand_minus_color,
        clip_on = True
    )
    ax.scatter(
        [ref_start, ref_end],
        [scf_start, scf_end] if strand == '+' else [scf_end, scf_start],
        s = 12,
        color = strand_plus_color if strand == '+' else strand_minus_color,
        clip_on = True
    )
# ax.set_xlim(0, ref_size)
ax.set_xlim(ref_size-1000000, ref_size)
# ax.set_ylim(0, scf_size)
ax.set_ylim(0, 1000000)
ax.set_xticks([])
ax.set_yticks([])
ax.set_xlabel(ref)
ax.set_ylabel(scf)
plt.show()
# %%
#####################
list_scf_order = table_paf_draw.sort_values(by = [9], ascending=False).drop_duplicates([0]).sort_values(by = [7])[0].to_list()
list_scf_order = sorted(list_scf_order)
# list_scf_order = ["ptg000023l", "ptg000063l", "ptg000005l", "ptg000011l"]
list_scf_order = ["ptg000029l", "ptg000047l", "ptg000005l", "ptg000043l"]
ref = "chr13"

# list_scf_order = list(filter(lambda name: dict_scf_to_scf_length[name] > 100000, list_scf_order))

#%%
from matplotlib import pyplot as plt

equal_scf_ratio = True

plt.rcParams["font.size"] = 20
fig, axes = plt.subplots(
    len(list_scf_order), 1,
    figsize=(8,20), 
    # figsize=(10,10), 
    gridspec_kw={
        'height_ratios': [1] * len(list_scf_order) if equal_scf_ratio else list(map(dict_scf_to_scf_length.__getitem__, list_scf_order))}
    # gridspec_kw={
    #     'height_ratios': list(map(dict_ref_to_ref_length.__getitem__, list_ref_order))[::-1],
    #     "width_ratios":list(map(dict_ref_to_ref_length.__getitem__, list_ref_order))}
)
plt.subplots_adjust(wspace = 0, hspace = 0)

strand_plus_color = "salmon"
strand_minus_color = "cornflowerblue"

for ind_scf, scf in enumerate(list_scf_order):
    axes[ind_scf].set_ylabel(scf, rotation = 0, ha = "right", va = 'center')
    if ind_scf == len(list_scf_order)-1:
        axes[ind_scf].set_xlabel(ref, rotation = 90)
    axes[ind_scf].set_xticks([])
    axes[ind_scf].set_yticks([])

    table_paf_draw_of_this_comb = table_paf_draw[
        np.logical_and(
            table_paf_draw[col_ref] == ref,
            table_paf_draw[col_qry] == scf
        )
    ]
    for _, row in table_paf_draw_of_this_comb.iterrows():
        scf_start = row[col_qry_start]
        scf_end = row[col_qry_end]
        scf_size = row[col_qry_len]
        ref_start = row[col_ref_start]
        ref_end = row[col_ref_end]
        ref_size = row[col_ref_len]
        strand = row[col_strand]
        assert strand in {'+', '-'}, row

        axes[ind_scf].plot(
            [ref_start, ref_end],
            [scf_start, scf_end] if strand == '+' else [scf_end, scf_start],
            linewidth = 4,
            color = strand_plus_color if strand == '+' else strand_minus_color,
            clip_on = True
        )
        axes[ind_scf].scatter(
            [ref_start, ref_end],
            [scf_start, scf_end] if strand == '+' else [scf_end, scf_start],
            s = 15,
            color = strand_plus_color if strand == '+' else strand_minus_color,
            clip_on = True
        )
        axes[ind_scf].set_xlim(0, ref_size)
        axes[ind_scf].set_xlim(0, 3000000)
        axes[ind_scf].set_ylim(0, scf_size)
plt.show()

# %%
###########################

from matplotlib import pyplot as plt

strand_plus_color = "salmon"
strand_minus_color = "cornflowerblue"

list_chr = list(map(lambda val: f"chr{val}", range(1, 23)))

plt.rcParams["font.size"] = 20
#%%
fig, axes = plt.subplots(
    4, 6,
    figsize=(24,16)
)
ax = axes.flatten()
for ind_chr, chrname in enumerate(list_chr):
    table_paf_draw_of_this_comb = table_paf_draw[
        np.logical_and(
            table_paf_draw[col_ref] == chrname,
            table_paf_draw[col_qry] == chrname
        )
    ]
    for _, row in table_paf_draw_of_this_comb.iterrows():
        scf_start = row[col_qry_start]
        scf_end = row[col_qry_end]
        scf_size = row[col_qry_len]
        ref_start = row[col_ref_start]
        ref_end = row[col_ref_end]
        ref_size = row[col_ref_len]
        strand = row[col_strand]
        assert strand in {'+', '-'}, row

        ax[ind_chr].plot(
            [ref_start, ref_end],
            [scf_start, scf_end] if strand == '+' else [scf_end, scf_start],
            linewidth = 3,
            color = strand_plus_color if strand == '+' else strand_minus_color,
            clip_on = True
        )
        ax[ind_chr].scatter(
            [ref_start, ref_end],
            [scf_start, scf_end] if strand == '+' else [scf_end, scf_start],
            s = 12,
            color = strand_plus_color if strand == '+' else strand_minus_color,
            clip_on = True
        )
        ax[ind_chr].set_xlim(0, ref_size)
        ax[ind_chr].set_ylim(0, scf_size)

    ax[ind_chr].set_xticks([])
    ax[ind_chr].set_yticks([])
    # ax[ind_chr].set_xlabel(chrname)
    # ax[ind_chr].set_ylabel(chrname)
    ax[ind_chr].set_title(chrname)
plt.show()
#%%
plt.rcParams["font.size"] = 15
plt.figure(figsize = (8, 3))
plt.bar(
    np.arange(len(list_chr))-0.2,
    list(map(lambda name: dict_ref_to_ref_length.get(name, 0), list_chr)),
    width = 0.4,
    label = "CHM13"
)
plt.bar(
    np.arange(len(list_chr))+0.2,
    list(map(lambda name: dict_scf_to_scf_length.get(name, 0), list_chr)),
    width = 0.4,
    label = "KOREF (Maternal)"
)
plt.xticks(np.arange(len(list_chr)), list_chr, rotation = 90)
plt.legend(loc = "upper center", bbox_to_anchor = (0.5, -0.4))
plt.xlim(-0.5, len(list_chr)-0.5)
plt.show()
# %%
