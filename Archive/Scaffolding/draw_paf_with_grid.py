#%%
import pandas as pd
import os
import numpy as np
from collections import Counter

hap = 1
# list_ref = ["chr13.Front", "chr13", "chr13.End", "chr15.Front", "chr15", "chr15.End"]
list_ref = ["chr13.Front1Mb", "chr13.End1Mb"]
list_qry = ["h1tg000041l","h1tg000044l","h1tg000060l"]
# list_qry = ["ntLink_3.Right.0", "h2tg000033l"]

path_paf_format = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling.RagTag.default_option.CHM13/DraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.telo_single.hap{hap}.RagTag.default_option.minimap2_latest.CHM13/OneToOneMapping/Alignment/Align.Between_Unmatched.Ref_CHM13.{{ref}}.Qry_{{qry}}.paf"

match_decision_ratio = 0.7
strand_decision_Ratio = 0.9

paf_filt_len = 0
paf_filt_quality = 0

#%%
list_table_paf = list()
for ref in list_ref:
    for qry in list_qry:
        table_paf_comb = pd.read_csv(path_paf_format.format(ref = ref, qry = qry), sep = '\t', names = list(range(24)))
        list_table_paf.append(table_paf_comb)

table_paf = pd.concat(list_table_paf)
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
# table_paf = table_paf[table_paf[col_match_len] >= paf_filt_len]
# table_paf = table_paf[table_paf[col_quality] >= paf_filt_quality]

# %%
dict_qry_to_ref_match = dict()

for qryname in list_qry:
    table_paf_qry = table_paf[table_paf[col_qry] == qryname]
    dict_qry_to_ref_match[qryname] = dict()
    for refname in list_ref:
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
table_paf_draw = table_paf[table_paf[col_match_len] >= paf_filt_len]
table_paf_draw = table_paf_draw[table_paf_draw[col_quality] >= paf_filt_quality]
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
list_scf_order = list_qry
list_ref_order = list_ref
#%%
from matplotlib import pyplot as plt

strand_plus_color = "salmon"
strand_minus_color = "cornflowerblue"

equal_scf_ratio = True

plt.rcParams["font.size"] = 25
fig, axes = plt.subplots(
    len(list_scf_order),len(list_ref_order),
    # figsize=(15,15), 
    figsize=(10,10), 
    gridspec_kw={
        'height_ratios': [1] * len(list_scf_order) if equal_scf_ratio else list(map(dict_scf_to_scf_length.__getitem__, list_scf_order)),
        "width_ratios":list(map(dict_ref_to_ref_length.__getitem__, list_ref_order))}
)
plt.subplots_adjust(wspace = 0.1, hspace = 0.1)

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
