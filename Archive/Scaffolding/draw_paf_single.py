#%%
import pandas as pd
import os
import numpy as np
from collections import Counter

path_paf = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/Local_Assembly/Hap2/chr15_98753195/Align/Align.REF_CHM13.QRY_assembly.paf"

match_decision_ratio = 0.7
strand_decision_Ratio = 0.9

paf_filt_len = 0
paf_filt_quality = 0

#%%
table_paf = pd.read_csv(path_paf, sep = '\t', names = list(range(24)))

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

table_paf = table_paf[table_paf[col_match_len] >= paf_filt_len]
table_paf = table_paf[table_paf[col_quality] >= paf_filt_quality]


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

table_paf_draw = table_paf[table_paf[col_match_len] >= paf_filt_len]
table_paf_draw = table_paf_draw[table_paf_draw[col_quality] >= paf_filt_quality]

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

list_scf_order = sorted(dict_scf_to_scf_length.keys(), key = dict_scf_to_scf_length.__getitem__)
list_ref_order = sorted(dict_ref_to_ref_length.keys(), key = dict_ref_to_ref_length.__getitem__)[::-1]
#%%
from matplotlib import pyplot as plt

strand_plus_color = "salmon"
strand_minus_color = "cornflowerblue"

###### Draw Single
ref = list_ref_order[0]
scf = list_scf_order[0]

ref_size = dict_ref_to_ref_length[ref]
scf_size = dict_scf_to_scf_length[scf]

show_size = min(ref_size, scf_size)
# show_size = show_size * 1.3

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

    ax.set_xlim(max(0, ref_size - show_size), ref_size)
    ax.set_ylim(0, min(show_size, scf_size))

ax.set_xticks([])
ax.set_yticks([])
ax.set_xlabel(ref)
ax.set_ylabel(scf)
plt.show()
# %%
