#%%
import os,sys,re,time
from copy import deepcopy
from joblib import Parallel, delayed
from collections import Counter
from pathlib import Path

import pandas as pd
import numpy as np

sys.path.append(str(Path(__file__).parents[2] / "common"))
from fasta_parser import read_fasta_as_dict, save_specific_seqnames_from_fasta_dict


def run(
    path_fasta,
    path_qv_perseq,
    path_hamming_perseq,
    path_qv_of_issued_region,
    path_lowqv_info,
    path_output,
    col_hamming_same,
    hamming_thresh = 0.05
):
    os.makedirs(os.path.dirname(path_output), exist_ok=True)

    dict_fasta = read_fasta_as_dict(path_fasta)
    table_qv = pd.read_csv(path_qv_perseq, names = ["name", "wrong_cnt", "total_cnt", "QV", "Err"], sep = '\t')
    table_qv_original = pd.read_csv(path_qv_of_issued_region, names = ["name", "wrong_cnt", "total_cnt", "QV", "Err"], sep = '\t')
    table_hamming = pd.read_csv(path_hamming_perseq, sep = '\t')
    table_lowqv_info = pd.read_csv(path_lowqv_info, sep = '\t')
    
    table_lowqv_info["regionname"] = table_lowqv_info.apply(lambda row: f"::{row['chr']}:{row['range_start']}-{row['range_end']}", axis = 1)
    table_lowqv_info["issueregionname"] = table_lowqv_info.apply(lambda row: f"::{row['chr']}:{row['lowqv_start']}-{row['lowqv_end']}", axis = 1)

    cols_hapmer = list(table_hamming.columns)[2:4]
    col_hamming_another = set(cols_hapmer) - {col_hamming_same}
    assert len(col_hamming_another) == 1
    col_hamming_another = list(col_hamming_another)[0]

    table_hamming["hamming_error"] = table_hamming[col_hamming_another] / (table_hamming[col_hamming_same] + table_hamming[col_hamming_another])
    table_hamming["hamming_error"] = table_hamming["hamming_error"].fillna(0)
    table_hamming_filt = table_hamming[table_hamming["hamming_error"] <= hamming_thresh]

    dict_issued_to_regionname = dict(zip(table_lowqv_info["issueregionname"], table_lowqv_info["regionname"]))
    table_qv_original["original_regionname"] = table_qv_original["name"].apply(dict_issued_to_regionname.__getitem__)
    dict_regionname_to_originalqv = dict(zip(table_qv_original["original_regionname"], table_qv_original["QV"]))
    
    table_qv["regionname"] = table_qv["name"].apply(lambda val: val.split("__")[0])
    table_qv["original_QV"] = table_qv["regionname"].apply(dict_regionname_to_originalqv.__getitem__)
    table_qv_filt = table_qv[table_qv["original_QV"] < table_qv["QV"]]

    contigids_get = set(table_qv_filt["name"]) & set(table_hamming_filt["Contig"])

    save_specific_seqnames_from_fasta_dict(dict_fasta, list(contigids_get), path_output)

#%%
if __name__ == "__main__":
    list_ref_assemblies = ['Hifiasm_Binned_HiFi', 'Hifiasm_Binned_ALT_HiFi', 'Hifiasm_Binned_HiFi_with0', 'Hifiasm_Binned_ALT_HiFi_with0', 'Flye_Binned_ONT', 'Flye_Binned_ONT.with0', 'Flye_Binned_HiFi_with0', 'Flye_Binned_HiFi']

    dict_hap_to_hamming_same = {
        1:"PATERNAL_Trimmed_F15_T5.noNbase.merged.k31.hapmer",
        2:"MATERNAL_Trimmed_F15_T5.noNbase.merged.k31.hapmer"
    }
    n_round = 7
    afterround = 4
    
    hamming_thresh = 1

    for hap in [1, 2]:
        path_qv_of_issued_region = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Manual_Polishing/Round{n_round}_afterR{afterround}/Collect_LowQV_Region/Merqury_QV/LowQVRegion.KOREF_S1.Trio_Specific.UpperCase.Polish_R{n_round-1}.Man.hap{hap}.QV_per_window.k31.chunked.2000.qv_less60.extra_range250.issue.sequence.k31.perseq.qv"
        path_lowqv_info = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Manual_Polishing/Round{n_round}_afterR{afterround}/Collect_LowQV_Region/LowQVRegion.KOREF_S1.Trio_Specific.UpperCase.Polish_R{n_round-1}.Man.hap{hap}.QV_per_window.k31.chunked.2000.qv_less60.extra_range250.region.tsv"
        
        for asmname in list_ref_assemblies:
            for anchor in ["with", "without"]:
                path_qv_perseq = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Manual_Polishing/Round{n_round}_afterR{afterround}/Anchored_Sequences/Filtered_Sequences/Hap{hap}/Merqury_QV/Anchored_Sequence.{asmname}.KOREF_S1.Trio_Specific.UpperCase.Polish_R{n_round-1}.Man.hap{hap}.sequence.without_anchor.k31.perseq.qv"
                path_hamming_perseq = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Manual_Polishing/Round{n_round}_afterR{afterround}/Anchored_Sequences/Filtered_Sequences/Hap{hap}/Merqury_Hamming/Anchored_Sequence.{asmname}.KOREF_S1.Trio_Specific.UpperCase.Polish_R{n_round-1}.Man.hap{hap}.sequence.with_anchor.hapmers.count"
                path_fasta = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Manual_Polishing/Round{n_round}_afterR{afterround}/Anchored_Sequences/Filtered_Sequences/Hap{hap}/Anchored_Sequence.{asmname}.KOREF_S1.Trio_Specific.UpperCase.Polish_R{n_round-1}.Man.hap{hap}.sequence.{anchor}_anchor.fasta"

                path_output = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Manual_Polishing/Round{n_round}_afterR{afterround}/Anchored_Sequences/Filtered_Sequences/Hap{hap}/Anchored_Sequence.{asmname}.KOREF_S1.Trio_Specific.UpperCase.Polish_R{n_round-1}.Man.hap{hap}.sequence.{anchor}_anchor.qv_improved.hamming{hamming_thresh}.fasta"
                
                assert os.path.exists(path_qv_perseq)
                assert os.path.exists(path_hamming_perseq)
                assert os.path.exists(path_fasta)
                assert os.path.exists(path_qv_of_issued_region)
                assert os.path.exists(path_lowqv_info)
                
                run(
                    path_fasta,
                    path_qv_perseq,
                    path_hamming_perseq,
                    path_qv_of_issued_region,
                    path_lowqv_info,
                    path_output,
                    dict_hap_to_hamming_same[hap],
                    hamming_thresh = hamming_thresh
                )
# %%
