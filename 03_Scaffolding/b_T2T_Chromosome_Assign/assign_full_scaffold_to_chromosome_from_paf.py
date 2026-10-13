#%%
import pandas as pd
import os
import numpy as np
from collections import Counter

hap = 2
refname = "CHM13"

path_paf = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling/Adjusted/Telomere_Detection/Fasta_Split_by_Telomere/Minimap2_Align/Align.KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.asm.dip.hap{hap}.p_ctg.no_telo.telo_single_masked.fa.k32.w100.z1000.ntLink.scaffolds.adjusted_telomere.telo_both.Ref_{refname}.paf"
path_save_assign = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling/Adjusted/Telomere_Detection/Fasta_Split_by_Telomere/Chromosome_Assign/Chromosome_Assign.KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.asm.dip.hap{hap}.p_ctg.no_telo.telo_single_masked.fa.k32.w100.z1000.ntLink.scaffolds.adjusted_telomere.telo_both.Ref_{refname}.tsv"

# path_paf = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling/Adjusted/Telomere_Detection/Fasta_Split_by_Telomere/Minimap2_Align/Align.KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg.no_telo.telo_single_masked.fa.k32.w100.z1000.ntLink.scaffolds.adjusted_telomere.telo_both.Ref_{refname}.paf"
# path_save_assign = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling/Adjusted/Telomere_Detection/Fasta_Split_by_Telomere/Chromosome_Assign/Chromosome_Assign.KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg.no_telo.telo_single_masked.fa.k32.w100.z1000.ntLink.scaffolds.adjusted_telomere.telo_both.Ref_{refname}.tsv"

match_decision_ratio = 0.7
strand_decision_Ratio = 0.9

paf_filt_len = 10000
paf_filt_quality = 30

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
# table_paf = table_paf[table_paf[col_match_len] >= paf_filt_len]
# table_paf = table_paf[table_paf[col_quality] >= paf_filt_quality]

# %%
dict_qry_to_ref_match = dict()

for qryname in table_paf[col_qry].unique():
    table_paf_qry = table_paf[table_paf[col_qry] == qryname]
    dict_qry_to_ref_match[qryname] = dict()
    for refname in table_paf_qry[col_ref].unique():
        table_paf_qry_to_ref = table_paf_qry[table_paf_qry[col_ref] == refname]

        reflen = table_paf_qry_to_ref.loc[:, col_ref_len].to_list()[0]
        all_matched_len = table_paf_qry_to_ref.loc[:, col_match_len].sum()

        strand_plus_len = table_paf_qry_to_ref[table_paf_qry_to_ref[col_strand] == '+'].loc[:, col_match_len].sum()
        strand_minus_len = table_paf_qry_to_ref[table_paf_qry_to_ref[col_strand] == '-'].loc[:, col_match_len].sum()

        if strand_plus_len > strand_minus_len:
            strand = '+'
            strand_ratio = strand_plus_len / (strand_plus_len + strand_minus_len)
        else:
            strand = '-'
            strand_ratio = strand_minus_len / (strand_plus_len + strand_minus_len)
        dict_qry_to_ref_match[qryname][refname] = {
            "matched_ratio" : all_matched_len/reflen,
            "strand" : strand,
            "strand_ratio" : strand_ratio
        }
# %%
dict_qry_to_ref_decision = dict()

for qryname in dict_qry_to_ref_match.keys():
    for refname in dict_qry_to_ref_match[qryname].keys():
        if dict_qry_to_ref_match[qryname][refname]["matched_ratio"] >= match_decision_ratio:
            assert dict_qry_to_ref_decision.get(qryname) == None, qryname
            assert dict_qry_to_ref_match[qryname][refname]["strand_ratio"] >= strand_decision_Ratio, dict_qry_to_ref_match[qryname][refname]

            dict_qry_to_ref_decision[qryname] = {
                "Matched_Reference" : refname,
                "Matched_Ratio" : dict_qry_to_ref_match[qryname][refname]["matched_ratio"],
                "Strand" : dict_qry_to_ref_match[qryname][refname]["strand"]
            }
# %%
table_qry_decision = pd.DataFrame(dict_qry_to_ref_decision).T
# %%
table_qry_decision_rename = table_qry_decision.reset_index(drop = False).rename(columns = {"index" : "Scaffold"})
# %%
os.makedirs(os.path.dirname(path_save_assign), exist_ok = True)
table_qry_decision_rename.to_csv(path_save_assign, sep = '\t', index = False)
# %%
#%%
table_qry_decision_rename_draw = table_qry_decision_rename.copy()
table_qry_decision_rename_draw["Scaffold_Length"] = table_qry_decision_rename_draw["Scaffold"].apply(lambda name: table_paf[table_paf[col_qry] == name].loc[:, col_qry_len].to_list()[0])
table_qry_decision_rename_draw["Matched_Reference_Length"] = table_qry_decision_rename_draw["Matched_Reference"].apply(lambda name: table_paf[table_paf[col_ref] == name].loc[:, col_ref_len].to_list()[0])

#%%
dict_chrname_to_rank = dict(zip(list(map(lambda val: f"chr{val}", range(1, 23)))+['chrX', 'chrY'], range(24)))
# %%
table_qry_decision_rename_draw["Ref_Rank"] = table_qry_decision_rename_draw["Matched_Reference"].apply(dict_chrname_to_rank.__getitem__)
table_qry_decision_rename_draw_sort = table_qry_decision_rename_draw.sort_values(by = "Ref_Rank", ascending=False)
# %%
from matplotlib import pyplot as plt

plt.rcParams["font.size"] = 14
plt.figure(figsize = (2, 2))

plt.barh(table_qry_decision_rename_draw_sort["Matched_Reference"], table_qry_decision_rename_draw_sort["Matched_Ratio"], color = table_qry_decision_rename_draw_sort["Strand"].apply({'+':"salmon", "-":"cornflowerblue"}.__getitem__))
plt.xlabel("Matched Ratio")
plt.ylim(-0.5, len(table_qry_decision_rename_draw_sort["Matched_Reference"])-0.5)
plt.axvline(0.8, linewidth = 2.5, color = "firebrick", linestyle = "--")
plt.show()
# %%
bar_width = 0.3

plt.rcParams["font.size"] = 14
plt.figure(figsize = (8, 4))

xtick_ref = np.arange(table_qry_decision_rename_draw_sort.shape[0])-bar_width/2
plt.bar(xtick_ref, table_qry_decision_rename_draw_sort["Matched_Reference_Length"][::-1], color = "gray", label = refname, width = bar_width)
xtick_scf = np.arange(table_qry_decision_rename_draw_sort.shape[0])+bar_width/2
plt.bar(xtick_scf, table_qry_decision_rename_draw_sort["Scaffold_Length"][::-1], color = "firebrick", label = "Scaffold", width = bar_width)

plt.xticks(range(table_qry_decision_rename_draw_sort.shape[0]), table_qry_decision_rename_draw_sort["Matched_Reference"][::-1], rotation = 90)
plt.legend(loc = "center left", bbox_to_anchor = (1.01, 0.5))
plt.ylabel("Length (bp)")
plt.xlabel("Chromosome")
plt.show()
# %%
table_paf_draw = table_paf[table_paf[col_match_len] >= paf_filt_len]
table_paf_draw = table_paf_draw[table_paf_draw[col_quality] >= paf_filt_quality]
# %%
dict_scf_to_scf_length = dict(zip(table_qry_decision_rename_draw_sort["Scaffold"], table_qry_decision_rename_draw_sort["Scaffold_Length"]))
dict_ref_to_ref_length = dict(zip(table_qry_decision_rename_draw_sort["Matched_Reference"], table_qry_decision_rename_draw_sort["Matched_Reference_Length"]))
# %%
list_scf_order = table_qry_decision_rename_draw_sort["Scaffold"]
list_ref_order = table_qry_decision_rename_draw_sort["Matched_Reference"][::-1]
#%%
strand_plus_color = "salmon"
strand_minus_color = "cornflowerblue"

plt.rcParams["font.size"] = 25
fig, axes = plt.subplots(
    len(list_scf_order),len(list_ref_order),
    figsize=(15,15), 
    gridspec_kw={
        'height_ratios': list(map(dict_scf_to_scf_length.__getitem__, list_scf_order)),
        "width_ratios":list(map(dict_ref_to_ref_length.__getitem__, list_ref_order))}
)
plt.subplots_adjust(wspace = 0, hspace = 0)

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
                color = strand_plus_color if strand == '+' else strand_minus_color
            )
            axes[ind_scf, ind_ref].set_xlim(0, ref_size)
            axes[ind_scf, ind_ref].set_ylim(0, scf_size)
plt.show()
# %%
