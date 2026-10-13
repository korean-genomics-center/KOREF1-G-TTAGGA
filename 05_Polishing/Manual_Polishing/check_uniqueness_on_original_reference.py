#%%
import os,sys
from copy import deepcopy
from joblib import Parallel, delayed
from collections import Counter
from pathlib import Path

import pandas as pd
import numpy as np

sys.path.append(str(Path(__file__).parents[2] / "common"))
from fasta_parser import read_fasta_as_dict, save_without_specific_seqnames_from_fasta_dict



def run(
    path_fasta,
    path_align,
    path_output,
    min_idy,
    min_align,
    max_hang,
    min_align2hang,
    min_rangein
):
    os.makedirs(os.path.dirname(path_output), exist_ok=True)
    table_paf = read_paf(path_align)
    table_paf_filt = filter_paf_for_uniqueness(table_paf, min_idy, min_align, max_hang, min_align2hang)

    list_seqnames_nonunique = get_seqnames_nonunique_on_desired(table_paf_filt, min_rangein)
    
    dict_fasta = read_fasta_as_dict(path_fasta)
    save_without_specific_seqnames_from_fasta_dict(
        dict_fasta,
        list_seqnames_nonunique,
        path_output
    )

def read_paf(path_paf):
    table_paf = pd.read_csv(path_paf, sep = '\t', names = list(range(1, 13)), usecols = list(range(1, 13)), low_memory = False)
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

def filter_paf_for_uniqueness(table_paf, min_idy, min_align, max_hang, min_align2hang):
    table_paf["identity"] = table_paf["match"] / table_paf["align"]

    table_paf_filt = table_paf[table_paf["identity"] > min_idy]
    table_paf_filt = table_paf_filt[table_paf_filt["align"] > min_align]

    table_paf_filt["wronghang"] = table_paf_filt.apply(get_wrong_hang_length, axis = 1)
    table_paf_filt["aligned2hang"] = table_paf_filt["align"] / (table_paf_filt["align"]+table_paf_filt["wronghang"])

    table_paf_filt = table_paf_filt[table_paf_filt["wronghang"] < max_hang]
    table_paf_filt = table_paf_filt[table_paf_filt["aligned2hang"] > min_align2hang]

    return table_paf_filt

def get_wrong_hang_length(dict_row):
    is_reverse_complement = dict_row["strand"] == '-'

    ref_left_hang = dict_row["ref_alignstart"]
    ref_right_hang = dict_row["reflen"] - dict_row["ref_alignend"]

    if is_reverse_complement:
        qry_left_hang = dict_row["qrylen"] - dict_row["qry_alignend"]
        qry_right_hang = dict_row["qry_alignstart"]
    else:
        qry_left_hang = dict_row["qry_alignstart"]
        qry_right_hang = dict_row["qrylen"] - dict_row["qry_alignend"]
    
    left_hang = min(ref_left_hang, qry_left_hang)
    right_hang = min(ref_right_hang, qry_right_hang)

    return left_hang + right_hang

def get_seqnames_nonunique_on_desired(table_paf_filt, min_rangein):
    list_seqnames_filtout = list()

    count_seqnames = Counter(table_paf_filt["qryname"])
    list_seqnames_nonunique = list(filter(lambda val: count_seqnames[val] >= 2, count_seqnames.keys()))
    list_seqnames_filtout.extend(list_seqnames_nonunique)

    list_seqnames_unique = list(filter(lambda val: count_seqnames[val] == 1, count_seqnames.keys()))
    for seqname_unique in list_seqnames_unique:
        regname = seqname_unique.split('__')[0].strip(':')
        chrname = regname.split(':')[0]
        start, end = regname.split(':')[1].split('-')
        start = int(start)
        end = int(end)

        table_paf_filt_of_seqname = table_paf_filt[table_paf_filt["qryname"] == seqname_unique]
        assert table_paf_filt_of_seqname.shape[0] == 1
        dict_align_info = table_paf_filt_of_seqname.iloc[0, :].to_dict()
        chr_aligned = dict_align_info["refname"]
        start_aligned = dict_align_info["ref_alignstart"]
        end_aligned = dict_align_info["ref_alignend"]

        if chrname != chr_aligned:
            list_seqnames_filtout.append(seqname_unique)
            continue

        n_range_overlap = len(set(range(start, end)) & set(range(start_aligned, end_aligned)))
        ratio_overlap = n_range_overlap / (end_aligned - start_aligned)
        if ratio_overlap < min_rangein:
            list_seqnames_filtout.append(seqname_unique)
    return list_seqnames_filtout



#%%
if __name__ == "__main__":
    min_idy = 0.95
    min_align = 10000
    max_hang = 1000
    min_align2hang = 0.95
    min_rangein = 0.95
    
    n_round = 7
    afterround = 4
    
    list_ref_assemblies = ['Hifiasm_Binned_HiFi', 'Hifiasm_Binned_ALT_HiFi', 'Hifiasm_Binned_HiFi_with0', 'Hifiasm_Binned_ALT_HiFi_with0', 'Flye_Binned_ONT', 'Flye_Binned_ONT.with0', 'Flye_Binned_HiFi_with0', 'Flye_Binned_HiFi']

    for hap in [1, 2]:
        for asmname in list_ref_assemblies:
            path_fasta = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Manual_Polishing/Round{n_round}_afterR{afterround}/Anchored_Sequences/Filtered_Sequences/Hap{hap}/Anchored_Sequence.{asmname}.KOREF_S1.Trio_Specific.UpperCase.Polish_R{n_round-1}.Man.hap{hap}.sequence.with_anchor.qv_improved.hamming1.fasta"

            path_align = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Manual_Polishing/Round{n_round}_afterR{afterround}/Anchored_Sequences/Filtered_Sequences/Hap{hap}/Minimap2_Align/Align.REF_PolishedASM_KOREF_S1_Polish_R{n_round-1}.Man.QRY_anchoredseq_from_{asmname}.with_anchor.qv_improved.hamming1.maphifiPw5s100cU101000000.filt_match1000_idy0.2.paf"

            path_output = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Manual_Polishing/Round{n_round}_afterR{afterround}/Anchored_Sequences/Filtered_Sequences.Unseen_Unique/Anchored_Sequence.{asmname}.KOREF_S1.Trio_Specific.UpperCase.Polish_R{n_round-1}.Man.hap{hap}.sequence.with_anchor.qv_improved.hamming1.unseen_unique.fasta"
            
            assert os.path.exists(path_fasta)
            assert os.path.exists(path_align)

            run(
                path_fasta,
                path_align,
                path_output,
                min_idy,
                min_align,
                max_hang,
                min_align2hang,
                min_rangein
            )
# %%
