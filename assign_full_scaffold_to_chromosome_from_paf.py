#%%
import pandas as pd
import os
import numpy as np
from collections import Counter

hap = 2
ref = "CHM13"

path_paf = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/Telomere_Detection/Fasta_Split_by_Telomere/Minimap2_Align/Align.KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.asm.dip.hap{hap}.p_ctg.telo_both.Ref_{ref}.paf"
path_save_assign = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/Telomere_Detection/Fasta_Split_by_Telomere/Chromosome_Assign/Chromosome_Assign.KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.asm.dip.hap{hap}.p_ctg.telo_both.Ref_{ref}.tsv"

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
dict_chrname_to_rank = dict(zip(list(map(lambda val: f"chr{val}", range(1, 23)))+['chrX', 'chrY'], range(24)))
# %%
table_qry_decision_rename["Ref_Rank"] = table_qry_decision_rename["Matched_Reference"].apply(dict_chrname_to_rank.__getitem__)
table_qry_decision_rename_sort = table_qry_decision_rename.sort_values(by = "Ref_Rank", ascending=False)
# %%
from matplotlib import pyplot as plt

plt.rcParams["font.size"] = 14
plt.figure(figsize = (2, 6))

plt.barh(table_qry_decision_rename_sort["Matched_Reference"], table_qry_decision_rename_sort["Matched_Ratio"], color = table_qry_decision_rename_sort["Strand"].apply({'+':"salmon", "-":"cornflowerblue"}.__getitem__))
plt.xlabel("Matched Ratio")
plt.ylim(-0.5, len(table_qry_decision_rename_sort["Matched_Reference"])-0.5)
plt.axvline(0.8, linewidth = 2.5, color = "firebrick", linestyle = "--")
plt.show()
# %%
