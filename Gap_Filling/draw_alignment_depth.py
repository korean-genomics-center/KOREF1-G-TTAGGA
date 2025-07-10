#%%
import os,sys
from copy import deepcopy

from joblib import Parallel, delayed
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt


def run(
    path_scaff_name,
    path_scaff_info,
    path_paf
):
    table_scaffname = pd.read_csv(path_scaff_name, sep = '\t', names = ["scaffname", "originalname"])
    dict_original2scaffname = dict(zip(table_scaffname["originalname"], table_scaffname["scaffname"].astype(str)))
    dict_scaff_info = read_scaff_info(path_scaff_info)

    table_paf = read_paf(path_paf)
    table_paf_highqual = table_paf[table_paf["align"] > 500]
    table_paf_highqual["identity"] = table_paf_highqual["match"] / table_paf_highqual["align"]
    table_paf_highqual = table_paf_highqual[table_paf_highqual["identity"] > 0.8]
    dict_ref_to_align_depth = get_alignment_depth_of_all_ref(table_paf_highqual)

def read_scaff_info(path_scaff_info):
    dict_scaff_info = dict()
    
    with open(path_scaff_info, 'r') as fr:
        for line in fr:
            if line.startswith('>'): continue
            values = line.strip('\n').split('\t')
            
            contigname = values[0]
            scaffname = values[6]
            contiglen = int(values[3])
            contigind = int(values[5])

            if dict_scaff_info.get(scaffname) == None:
                dict_scaff_info[scaffname] = dict()
            assert dict_scaff_info[scaffname].get(contigind) == None
            dict_scaff_info[scaffname][contigind] = {
                "name" : contigname,
                "length" : contiglen
            }
    return dict_scaff_info

def read_paf(path_paf):
    table_paf = pd.read_csv(path_paf, sep = '\t', names = list(range(1, 13)), usecols = list(range(1, 13)))
    table_paf = table_paf.rename(
        columns = {
            1:"qryname",
            2:"qrylen",
            3:"qry_alignstart",
            4:"qry_alignend",
            5:"strand",
            6:"refname",
            7:"reflen",
            8:"ref_alignstart",
            9:"ref_alignend",
            10:"match",
            11:"align",
            12:"quality"
        }
    )
    return table_paf

def get_alignment_depth_of_all_ref(table_paf):
    dict_ref_to_length = get_ref_to_length_from_paf(table_paf)
    list_ref = sorted(dict_ref_to_length.keys())

    table_paf_indexed = table_paf.set_index("refname", drop = False)

    with Parallel(n_jobs=len(dict_ref_to_length)) as parallel:
        list_depth_per_ref = parallel(
            delayed(get_alignment_cumulative_count_of_specific_ref)(
                table_paf_indexed.loc[refname, :], dict_ref_to_length[refname]
            ) for refname in list_ref
        )
    dict_ref_to_align_depth = dict(zip(list_ref, list_depth_per_ref))

    return dict_ref_to_align_depth

def get_ref_to_length_from_paf(table_paf):
    table_paf_ref_uniq = table_paf.drop_duplicates(subset=["refname"])

    dict_ref_to_length = dict(zip(table_paf_ref_uniq["refname"], table_paf_ref_uniq["reflen"]))
    return dict_ref_to_length

def get_alignment_cumulative_count_of_specific_ref(table_paf, ref_len):
    arr_values = np.zeros(ref_len+1)

    for _, row in table_paf.iterrows():
        ref_start = row["ref_alignstart"]
        ref_end = row["ref_alignend"]

        arr_values[np.arange(ref_start, ref_end+1)] += 1
    return arr_values

def init_count_dictionary_from_length(dict_name_to_length):
    dict_count = dict()
    for name, length in dict_name_to_length.items():
        dict_count[name] = dict(zip(
            range(1, length+1),
            [0] * length
        ))
    return dict_count



# %%
if __name__ == "__main__":
    path_scaff_name = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon/Round13.TGS_Gapcloser_Modified_V2_re.Unpl_Scf.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.Manual_Scaff.HiFi_trio.Hifiasm_Loc.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.ONT_trio.ONT_Amb.Manual_Scaff/Prealignment_for_Candidate_Filtering/Hap1/CuratedDraftASM.KOREF_S1.Trio_Specific.hap1.GapFilling.Unplaced_Scf.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.Manual_Scaff.HiFi_trio.Hifiasm_Loc.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.ONT_trio.ONT_Amb.Manual_Scaff.name_map"
    path_scaff_info = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon/Round13.TGS_Gapcloser_Modified_V2_re.Unpl_Scf.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.Manual_Scaff.HiFi_trio.Hifiasm_Loc.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.ONT_trio.ONT_Amb.Manual_Scaff/Prealignment_for_Candidate_Filtering/Hap1/CuratedDraftASM.KOREF_S1.Trio_Specific.hap1.GapFilling.Unplaced_Scf.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.Manual_Scaff.HiFi_trio.Hifiasm_Loc.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.ONT_trio.ONT_Amb.Manual_Scaff.orignial_scaff_infos"
    path_paf = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon/Round13.TGS_Gapcloser_Modified_V2_re.Unpl_Scf.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.Manual_Scaff.HiFi_trio.Hifiasm_Loc.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.ONT_trio.ONT_Amb.Manual_Scaff/Prealignment_for_Candidate_Filtering/Hap1/Align/Align.REF_KOREFDraftASM_Hap1_Split.QRY_HiFi_all.maphifi.paf"
# %%
from matplotlib import pyplot as plt

def get_mean_window_value(list_values, window_size = 1000):
    list_windowed = list()
    start = 0
    while 1:
        values_part = list_values[start:start+window_size]
        list_windowed.append(np.mean(values_part))
        if (start+window_size) >= len(list_values):
            break
        start = start+window_size
    return list_windowed

window_size = 1000

for key in dict_ref_to_align_depth.keys():
    print(key)
    list_values = dict_ref_to_align_depth[key]
    list_windowed_values = get_mean_window_value(list_values, window_size)

    plt.rcParams["font.size"] = 15
    plt.figure(figsize = (20, 3))
    plt.scatter(
        np.arange(len(list_windowed_values))*window_size,
        list_windowed_values,
        s = 5,
        color = "gray"
    )
    plt.xlabel("Coordinate (bp)")
    plt.ylabel("Depth")
    xlim_pad = len(list_values)/50
    # plt.xlim(len(list_values)-1e6, len(list_values)+xlim_pad)
    plt.xlim(-xlim_pad, len(list_values)+xlim_pad)
    plt.ylim(top = np.quantile(list_values, 0.999))
    plt.axvline(0, linewidth = 1, linestyle = "--", color = "k")
    plt.axvline(len(list_values), linewidth = 1, linestyle = "--", color = "k")
    plt.show()
# %%
