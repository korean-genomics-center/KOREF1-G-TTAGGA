#%%
import os,sys
from copy import deepcopy
from joblib import Parallel, delayed
from collections import Counter
from pathlib import Path

import pandas as pd
import numpy as np

sys.path.append(str(Path(__file__).parents[2]))
from fasta_parser import read_fasta_as_dict, save_dict_fasta_into_file


def run(
    path_original_ref,
    path_lowqv_region,
    path_selected_fasta,
    path_selected_info,
    path_output
):
    os.makedirs(os.path.dirname(path_output), exist_ok=True)

    dict_original_ref = read_fasta_as_dict(path_original_ref)
    dict_selected_fasta = read_fasta_as_dict(path_selected_fasta)

    table_lowqv = pd.read_csv(path_lowqv_region, sep = '\t')
    table_selected_info = pd.read_csv(path_selected_info, sep = '\t')
    
    table_lowqv["regionname"] = table_lowqv.apply(lambda row: f"::{row['chr']}:{row['range_start']}-{row['range_end']}", axis = 1)
    table_lowqv_indexed = table_lowqv.set_index("regionname")
    table_selected_info_indexed = table_selected_info.set_index("region", drop = False)

    dict_replaced_ref = replace_original_reference_with_new_sequence(dict_original_ref, dict_selected_fasta, table_lowqv_indexed, table_selected_info_indexed)

    save_dict_fasta_into_file(dict_replaced_ref, path_output)

def replace_original_reference_with_new_sequence(dict_original_ref, dict_selected_fasta, table_lowqv_indexed, table_selected_info_indexed):
    dict_replaced_ref = deepcopy(dict_original_ref)

    table_lowqv_indexed_selected = table_lowqv_indexed.loc[table_selected_info_indexed["region"].to_list()]
    table_lowqv_indexed_selected_sort = table_lowqv_indexed_selected.sort_values(by = ["lowqv_end"], ascending = False)

    for regionname, row in table_lowqv_indexed_selected_sort.iterrows():
        chrname = row["chr"]
        lowqv_start = row["lowqv_start"]
        lowqv_end = row["lowqv_end"]

        dict_selected_thisone = table_selected_info_indexed.loc[regionname, :].to_dict()

        sequence_for_patching = dict_selected_fasta[dict_selected_thisone["seqname"]]
        assert len(sequence_for_patching) == (dict_selected_thisone["front_anchor_len"] + dict_selected_thisone["tail_anchor_len"] + dict_selected_thisone["replace_len"])

        sequence_for_patching_wo_anchor = sequence_for_patching[dict_selected_thisone["front_anchor_len"]:-dict_selected_thisone["tail_anchor_len"]]

        sequence_ref_front = dict_replaced_ref[chrname][:lowqv_start]
        sequence_ref_end = dict_replaced_ref[chrname][lowqv_end:]

        dict_replaced_ref[chrname] = sequence_ref_front+sequence_for_patching_wo_anchor+sequence_ref_end
    return dict_replaced_ref

#%%
if __name__ == "__main__":
    n_round = 7
    afterround = 4
    
    for hap in [1, 2]:

        path_original_ref = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Reference_Fasta/Manual_Polishing/After_Round{afterround}/Round{n_round-1}.Manual_Polished_Fasta/PolishedASM.KOREF_S1.Trio_Specific.UpperCase.Polish_R{n_round-1}.Man.hap{hap}.fasta"
        path_lowqv_region = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Manual_Polishing/Round{n_round}_afterR{afterround}/Collect_LowQV_Region/LowQVRegion.KOREF_S1.Trio_Specific.UpperCase.Polish_R{n_round-1}.Man.hap{hap}.QV_per_window.k31.chunked.2000.qv_less60.extra_range250.region.tsv"
        path_selected_fasta = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Manual_Polishing/Round{n_round}_afterR{afterround}/Anchored_Sequence_for_Polishing/Anchored_Sequence.KOREF_S1.Trio_Specific.UpperCase.Polish_R{n_round-1}.Man.hap{hap}.sequence.with_anchor.qv_improved.hamming1.unseen_unique.selected.fasta"
        path_selected_info = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Manual_Polishing/Round{n_round}_afterR{afterround}/Anchored_Sequence_for_Polishing/Anchored_Sequence.KOREF_S1.Trio_Specific.UpperCase.Polish_R{n_round-1}.Man.hap{hap}.sequence.with_anchor.qv_improved.hamming1.unseen_unique.selected.info.tsv"

        path_output = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Reference_Fasta/Manual_Polishing/After_Round{afterround}/Round{n_round}.Manual_Polished_Fasta/PolishedASM.KOREF_S1.Trio_Specific.UpperCase.Polish_R{n_round}.Man.hap{hap}.fasta"

        run(
            path_original_ref,
            path_lowqv_region,
            path_selected_fasta,
            path_selected_info,
            path_output
        )
# %%
