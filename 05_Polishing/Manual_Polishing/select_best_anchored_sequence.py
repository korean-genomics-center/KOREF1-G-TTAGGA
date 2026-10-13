#%%
import os,sys
from copy import deepcopy
from joblib import Parallel, delayed
from collections import Counter
from pathlib import Path

import pandas as pd
import numpy as np

sys.path.append(str(Path(__file__).parents[2] / "common"))
from fasta_parser import read_fasta_as_dict, save_dict_fasta_into_file


def run(
    path_fasta_format,
    path_anchor_info_format,
    path_qv_format,
    list_asmnames,
    path_output_prefix
):
    os.makedirs(os.path.dirname(path_output_prefix), exist_ok=True)
    dict_name_to_fasta = read_all_fasta_by_format(path_fasta_format, list_asmnames)
    dict_region_to_filtedanchor_info = read_and_filter_anchorinfo_by_fasta(dict_name_to_fasta, path_anchor_info_format, path_qv_format, list_asmnames)
    dict_region_to_selected_seq = select_best_sequence_from_table(dict_region_to_filtedanchor_info)
    dict_selected_fasta = get_fasta_of_selected_seq(dict_name_to_fasta, dict_region_to_selected_seq)

    table_info_all = pd.concat(dict_region_to_filtedanchor_info.values())
    table_info_selected = table_info_all[table_info_all["seqname"].isin(list(dict_selected_fasta.keys()))]

    path_save_fasta = f"{path_output_prefix}.fasta"
    path_save_info = f"{path_output_prefix}.info.tsv"

    save_dict_fasta_into_file(dict_selected_fasta, path_save_fasta)
    table_info_selected.to_csv(path_save_info, sep = '\t', index = False)

def read_all_fasta_by_format(path_fasta_format, list_names):
    dict_name_to_fasta = dict()

    for name in list_names:
        dict_name_to_fasta[name] = read_fasta_as_dict(path_fasta_format.format(asmname = name))
    return dict_name_to_fasta

def read_and_filter_anchorinfo_by_fasta(dict_name_to_fasta, path_anchor_info_format, path_qv_format, list_amsnames):
    dict_region_to_anchorinfo_list = dict()

    for name in list_amsnames:
        table_anchorinfo_of_name = pd.read_csv(path_anchor_info_format.format(asmname = name), sep = '\t')
        table_qv_of_name = pd.read_csv(path_qv_format.format(asmname = name), names = ["name", "wrong_cnt", "total_cnt", "QV", "Err"], sep = '\t')
        dict_seqname_to_qv = dict(zip(table_qv_of_name["name"], table_qv_of_name["QV"]))

        for seqname in dict_name_to_fasta[name].keys():
            regionname = seqname.split("__")[0]
            contigname = seqname.split("__")[-1]

            if dict_region_to_anchorinfo_list.get(regionname) == None:
                dict_region_to_anchorinfo_list[regionname] = list()

            table_anchorinfo_of_name_this_seqname = table_anchorinfo_of_name[np.logical_and(table_anchorinfo_of_name["region"] == regionname, table_anchorinfo_of_name["contig"] == contigname)]            
            assert table_anchorinfo_of_name_this_seqname.shape[0] == 1
            table_anchorinfo_of_name_this_seqname["asmname"] = name
            table_anchorinfo_of_name_this_seqname["seqname"] = seqname
            table_anchorinfo_of_name_this_seqname["QV"] = table_anchorinfo_of_name_this_seqname["seqname"].apply(dict_seqname_to_qv.__getitem__)

            dict_region_to_anchorinfo_list[regionname].append(table_anchorinfo_of_name_this_seqname.copy())
    
    dict_region_to_anchorinfo = dict()
    for regionname, list_info in dict_region_to_anchorinfo_list.items():
        dict_region_to_anchorinfo[regionname] = pd.concat(list_info)
    return dict_region_to_anchorinfo

def select_best_sequence_from_table(dict_region_to_filtedanchor_info):
    dict_region_to_selected = dict()

    for regionname, table_filteredanchor_info in dict_region_to_filtedanchor_info.items():
        best_qv = table_filteredanchor_info["QV"].max()
        table_filteredanchor_info_bestqv = table_filteredanchor_info[table_filteredanchor_info["QV"] == best_qv]
        
        table_filteredanchor_info_bestqv["mean_anchor_idy"] = (table_filteredanchor_info_bestqv["front_anchor_idy"] + table_filteredanchor_info_bestqv["tail_anchor_idy"]) / 2

        best_mean_idy = table_filteredanchor_info_bestqv["mean_anchor_idy"].max()
        table_filteredanchor_info_bestidy = table_filteredanchor_info_bestqv[table_filteredanchor_info_bestqv["mean_anchor_idy"] == best_mean_idy]

        table_filteredanchor_info_bestidy["total_anchor_len"] = table_filteredanchor_info_bestidy["front_anchor_len"] + table_filteredanchor_info_bestidy["tail_anchor_len"]
        table_filteredanchor_info_bestidy = table_filteredanchor_info_bestidy.sort_values(by = ["total_anchor_len"], ascending = False)

        row_selected = table_filteredanchor_info_bestidy.iloc[0, :].to_dict()

        dict_region_to_selected[regionname] = {
            "asmname" : row_selected["asmname"],
            "seqname" : row_selected["seqname"]
        }
    return dict_region_to_selected

def get_fasta_of_selected_seq(dict_name_to_fasta, dict_region_to_selected_seq):
    dict_selected_fasta = dict()

    for regionname in dict_region_to_selected_seq.keys():
        seqname = dict_region_to_selected_seq[regionname]["seqname"]
        dict_selected_fasta[seqname] = dict_name_to_fasta[dict_region_to_selected_seq[regionname]["asmname"]][seqname]
    return dict_selected_fasta

#%%
if __name__ == "__main__":
    n_round = 7
    afterround = 4
    
    list_asmnames = ['Hifiasm_Binned_HiFi', 'Hifiasm_Binned_ALT_HiFi', 'Hifiasm_Binned_HiFi_with0', 'Hifiasm_Binned_ALT_HiFi_with0', 'Flye_Binned_ONT', 'Flye_Binned_ONT.with0', 'Flye_Binned_HiFi_with0', 'Flye_Binned_HiFi']

    for hap in [1, 2]:
        path_fasta_format = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Manual_Polishing/Round{n_round}_afterR{afterround}/Anchored_Sequences/Filtered_Sequences.Unseen_Unique/Anchored_Sequence.{{asmname}}.KOREF_S1.Trio_Specific.UpperCase.Polish_R{n_round-1}.Man.hap{hap}.sequence.with_anchor.qv_improved.hamming1.unseen_unique.fasta"
        path_anchor_info_format = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Manual_Polishing/Round{n_round}_afterR{afterround}/Anchored_Sequences/Filtered_Sequences/Hap{hap}/Anchored_Sequence.{{asmname}}.KOREF_S1.Trio_Specific.UpperCase.Polish_R{n_round-1}.Man.hap{hap}.info.tsv"
        path_qv_format = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Manual_Polishing/Round{n_round}_afterR{afterround}/Anchored_Sequences/Filtered_Sequences/Hap{hap}/Merqury_QV/Anchored_Sequence.{{asmname}}.KOREF_S1.Trio_Specific.UpperCase.Polish_R{n_round-1}.Man.hap{hap}.sequence.without_anchor.k31.perseq.qv"

        path_output_prefix = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Manual_Polishing/Round{n_round}_afterR{afterround}/Anchored_Sequence_for_Polishing/Anchored_Sequence.KOREF_S1.Trio_Specific.UpperCase.Polish_R{n_round-1}.Man.hap{hap}.sequence.with_anchor.qv_improved.hamming1.unseen_unique.selected"
        
        run(path_fasta_format, path_anchor_info_format, path_qv_format, list_asmnames, path_output_prefix)

# %%
