#%%
import os

import pandas as pd

#%%
dict_metrics_format = {
    "Length" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Final_Assembly/FinalASM.KOREF_S1.Trio_Specific.hap{hap}.Reordered.fasta.fai",
    "QV" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Merqury_Evaluation/Merqury_QV/PolishedASM.KOREF_S1.Trio_Specific.UpperCase.Polish_R8.SV_SNV.hap{hap}.k21.perseq.qv",
    "GCI" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/GCI_Evaluation/GCI_Result/Hap{hap}/GCI_Result.hap{hap}.HiFi_ONT.Winnowmap2_Minimap2.triobin_{triobin}0a.gci",
    "CRAQ" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/CRAQ_Evaluation/CRAQ_Result/LongRead_Merged_All.ShortRead/Hap{hap}/runAQI_out/out_final.Report",
    "Flagger" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Flagger_Evaluation/Flagger/ONT_10kb_QV15/final_flagger_prediction.Hap{hap}.bed"
}
# %%
def read_length_from_fai(path_fai):
    dict_chrname_to_length = dict()
    with open(path_fai, 'r') as fr:
        for line in fr:
            values = line.strip('\n').split('\t')
            chrname, length = values[:2]
            dict_chrname_to_length[chrname] = int(length)
    return dict_chrname_to_length

def read_merqury_qv(path_qv_per_seq):
    dict_chrname_to_qv = dict()
    with open(path_qv_per_seq, 'r') as fr:
        for line in fr:
            values = line.strip('\n').split('\t')
            chrname = values[0]
            qv = values[3]
            if qv == "inf":
                qv = 100
            else:
                qv = float(qv)
            dict_chrname_to_qv[chrname] = qv
    return dict_chrname_to_qv

def read_gci_score(path_gci_result):
    dict_chrname_to_gci = dict()
    with open(path_gci_result, 'r') as fr:
        is_combined_score_line_started = False
        is_header_passed = False

        for line in fr:
            if line.strip('\n') == "HiFi + Nano:":
                is_combined_score_line_started = True
                continue
            if is_combined_score_line_started:
                if not is_header_passed:
                    is_header_passed = True
                    continue
                values = line.strip('\n').split('\t')
                if len(values) != 6:
                    continue
                else:
                    chrname = values[0]
                    gci_score = float(values[5])
                    dict_chrname_to_gci[chrname] = gci_score
    return dict_chrname_to_gci

def read_craq_aqi_score(path_craq_out):
    dict_chrname_to_aqi = dict()
    with open(path_craq_out, 'r') as fr:
        for line in fr:
            values = line.strip('\n').split('\t')
            if len(values) != 8: continue
            if values[0] == "#Chr" : continue
            chrname = values[0]
            aqi_score = float(values[7])
            dict_chrname_to_aqi[chrname] = aqi_score
    return dict_chrname_to_aqi

def read_flagger_result(path_flagger_bed, chrname_suffix, dict_chrname_to_length):
    dict_flagger_region = dict()
    with open(path_flagger_bed, 'r') as fr:
        for line in fr:
            list_values = line.strip('\n').split('\t')
            if len(list_values) < 3: continue

            chrname, start, end, flag = list_values[:4]
            chrname = chrname.replace(chrname_suffix, '')

            if dict_flagger_region.get(chrname) == None:
                dict_flagger_region[chrname] = 0
            
            if flag != "Hap":
                dict_flagger_region[chrname] += (int(end) - int(start))
    
    dict_flagger_ratio = dict()
    for chrname, chr_length in dict_chrname_to_length.items():
        dict_flagger_ratio[chrname] = 1 - dict_flagger_region[chrname] / chr_length
    return dict_flagger_ratio


# %%
dict_hap_to_trio = {
    1 : "p",
    2 : "m"
}

dict_metric_to_func = {
    "Length" : read_length_from_fai,
    "QV" : read_merqury_qv,
    "GCI" : read_gci_score,
    "CRAQ" : read_craq_aqi_score,
    "Flagger" : read_flagger_result
}

dict_hap_to_chrs = {
    1 : list(map(lambda val: f"chr{val}", range(1, 23))) + ["chrY"],
    2 : list(map(lambda val: f"chr{val}", range(1, 23))) + ["chrX"],
}

dict_hap_to_metrics = dict()
for hap in dict_hap_to_chrs.keys():
    dict_hap_to_metrics[hap] = dict()
    for metric, path_metric_format in dict_metrics_format.items():
        dict_hap_to_metrics[hap][metric] = dict()
        path_metric = path_metric_format.format(hap = hap, triobin = dict_hap_to_trio[hap])

        if metric == "Flagger":
            dict_metric_value = dict_metric_to_func[metric](path_metric, f"_Hap{hap}", dict_hap_to_metrics[hap]["Length"])
        else:
            dict_metric_value = dict_metric_to_func[metric](path_metric)
        
        for chrname in dict_hap_to_chrs[hap]:
            dict_hap_to_metrics[hap][metric][chrname] = dict_metric_value[chrname]
# %%
table_hap1_metric = pd.DataFrame(dict_hap_to_metrics[1])
table_hap2_metric = pd.DataFrame(dict_hap_to_metrics[2])
# %%
table_hap1_metric.columns = list(map(lambda val: f"{val}_Hap1", table_hap1_metric.columns))
table_hap2_metric.columns = list(map(lambda val: f"{val}_Hap2", table_hap2_metric.columns))
# %%
table_metric_total = pd.concat([table_hap1_metric, table_hap2_metric], axis = 1)
# %%
dict_vote = dict()
for chrname in list(map(lambda val: f"chr{val}", range(1, 23))) + ["chrY", "chrX"]:
    dict_vote[chrname] = dict()
    for metric in dict_metrics_format.keys():
        metric_hap1 = table_metric_total.loc[chrname, f"{metric}_Hap1"]
        metric_hap2 = table_metric_total.loc[chrname, f"{metric}_Hap2"]

        vote = None
        if pd.isna(metric_hap1):
            vote = "Maternal"
        elif pd.isna(metric_hap2):
            vote = "Paternal"
        elif metric_hap1 == metric_hap2:
            vote = "Equal"
        else:
            vote = "Paternal" if metric_hap1 > metric_hap2 else "Maternal"
        
        dict_vote[chrname][f"{metric}_Vote"] = vote
        

# %%
table_vote = pd.DataFrame(dict_vote).T
# %%
def vote_bestone(dict_values):
    from collections import Counter
    count_values = Counter(dict_values.values())
    if count_values["Maternal"] > count_values["Paternal"]:
        return "Maternal"
    elif count_values["Maternal"] < count_values["Paternal"]:
        return "Paternal"
    else:
        return dict_values["Length_Vote"]

table_vote["Combined"] = table_vote.apply(lambda row: vote_bestone(row.to_dict()), axis = 1)
# %%
table_metric_total_with_vote = pd.concat([table_metric_total, table_vote], axis = 1)
# %%
table_metric_total_with_vote.to_csv("/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Final_Assembly/Vote_Chromosomes_for_Combined.tsv", sep = '\t')
# %%
