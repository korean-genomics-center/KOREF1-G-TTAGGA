#%%
import os,sys,re,time
from copy import deepcopy
from joblib import Parallel, delayed
from collections import Counter
from pathlib import Path

import pandas as pd
import numpy as np

sys.path.append(str(Path(__file__).parents[2]))
from fasta_parser import read_fasta_as_dict, save_dict_fasta_into_file


def run(
    path_alignment_on_anchor,
    path_original_fasta,
    path_output_prefix,
    name
):
    os.makedirs(os.path.dirname(path_output_prefix), exist_ok=True)
    table_alignment_on_anchor = pd.read_csv(path_alignment_on_anchor, sep = '\t')
    dict_anchor_comb_per_qry = collect_all_anchor_combinations(table_alignment_on_anchor)
    dict_original_fasta = read_fasta_as_dict(path_original_fasta)

    dict_anchored_sequence = get_anchored_contig_sequence(dict_anchor_comb_per_qry, dict_original_fasta, name)

    path_output_seq_with_anchor = f"{path_output_prefix}.sequence.with_anchor.fasta"
    path_output_seq_without_anchor = f"{path_output_prefix}.sequence.without_anchor.fasta"
    path_output_seq_info = f"{path_output_prefix}.info.tsv"
    save_anchored_sequence_with_info(dict_anchored_sequence, path_output_seq_with_anchor, path_output_seq_without_anchor, path_output_seq_info)

def collect_all_anchor_combinations(table_alignment_on_anchor):
    list_query = table_alignment_on_anchor["qryname"].unique()

    dict_anchor_comb_per_qry = dict()
    for qry in list_query:
        table_alignment_on_anchor_single = table_alignment_on_anchor[table_alignment_on_anchor["qryname"] == qry]

        dict_anchors_on_single_qry = check_bothside_anchors_of_single_query(table_alignment_on_anchor_single)

        if len(dict_anchors_on_single_qry) > 0:
            dict_anchor_comb_per_qry[qry] = dict_anchors_on_single_qry
    return dict_anchor_comb_per_qry

def check_bothside_anchors_of_single_query(table_single_qry):
    dict_refalignment_of_anchors = dict()
    for refname in table_single_qry["refname"].unique():
        table_single_qry_single_ref = table_single_qry[table_single_qry["refname"] == refname]

        table_single_qry_single_ref_bothanchor = table_single_qry_single_ref[np.logical_and(table_single_qry_single_ref["front_anchor"], table_single_qry_single_ref["tail_anchor"])]
        ind_bothanchor = list(table_single_qry_single_ref_bothanchor.index)

        dict_both_anchor = None
        if table_single_qry_single_ref_bothanchor.shape[0] > 0:
            table_single_qry_single_ref_bothanchor["mean_identity"] = (table_single_qry_single_ref_bothanchor["front_identity"] + table_single_qry_single_ref_bothanchor["tail_identity"]) / 2
            
            dict_both_anchor = table_single_qry_single_ref_bothanchor.sort_values(by = ["mean_identity"], ascending = False).iloc[0, :].to_dict()

        table_single_qry_single_ref_singleanchor = table_single_qry_single_ref.loc[list(filter(lambda ind: ind not in ind_bothanchor, table_single_qry_single_ref.index))]

        table_single_qry_single_ref_frontanchor = table_single_qry_single_ref_singleanchor[table_single_qry_single_ref_singleanchor["front_anchor"]]
        table_single_qry_single_ref_tailanchor = table_single_qry_single_ref_singleanchor[table_single_qry_single_ref_singleanchor["tail_anchor"]]

        if table_single_qry_single_ref_frontanchor.shape[0] == 0 or table_single_qry_single_ref_tailanchor.shape[0] == 0:
            if dict_both_anchor != None:
                dict_refalignment_of_anchors[refname] = {
                    "front":dict_both_anchor,
                    "tail":dict_both_anchor
                }
            continue

        dict_front_anchor = table_single_qry_single_ref_frontanchor.sort_values(by = ["front_identity"], ascending = False).iloc[0, :].to_dict()
        dict_tail_anchor = table_single_qry_single_ref_tailanchor.sort_values(by = ["tail_identity"], ascending = False).iloc[0, :].to_dict()

        if dict_front_anchor["strand"] != dict_tail_anchor["strand"]:
            if dict_both_anchor != None:
                dict_refalignment_of_anchors[refname] = {
                    "front":dict_both_anchor,
                    "tail":dict_both_anchor
                }
            continue
        if dict_front_anchor["strand"] == '+':
            ref_frontanchor = dict_front_anchor["ref_alignstart"] + dict_front_anchor["front_anchor_alignlen"] - 1
            ref_tailanchor = dict_tail_anchor["ref_alignend"] - dict_tail_anchor["tail_anchor_alignlen"] + 1
            if ref_frontanchor > ref_tailanchor:
                if dict_both_anchor != None:
                    dict_refalignment_of_anchors[refname] = {
                        "front":dict_both_anchor,
                        "tail":dict_both_anchor
                    }
                continue
        else:
            ref_frontanchor = dict_front_anchor["ref_alignend"] - dict_front_anchor["front_anchor_alignlen"] + 1
            ref_tailanchor = dict_tail_anchor["ref_alignstart"] + dict_tail_anchor["tail_anchor_alignlen"] - 1
            if ref_frontanchor < ref_tailanchor:
                if dict_both_anchor != None:
                    dict_refalignment_of_anchors[refname] = {
                        "front":dict_both_anchor,
                        "tail":dict_both_anchor
                    }
                continue

        if dict_both_anchor != None and (dict_front_anchor["front_identity"] + dict_tail_anchor["tail_identity"]) / 2 < dict_both_anchor["mean_identity"]:
            dict_refalignment_of_anchors[refname] = {
                "front":dict_both_anchor,
                "tail":dict_both_anchor
            }
        else:
            dict_refalignment_of_anchors[refname] = {
                "front":dict_front_anchor,
                "tail":dict_tail_anchor
            }
    return dict_refalignment_of_anchors

def get_anchored_contig_sequence(dict_anchor_comb_per_qry, dict_original_fasta, name):
    dict_sequence_info = dict()

    for regname in dict_anchor_comb_per_qry.keys():
        for contigname, dict_anchor_align_single in dict_anchor_comb_per_qry[regname].items():
            strand = dict_anchor_align_single["front"]["strand"]
            if strand == '+':
                contig_start = dict_anchor_align_single["front"]["ref_alignstart"]
                contig_end = dict_anchor_align_single["tail"]["ref_alignend"]

                anchored_sequence = dict_original_fasta[contigname][contig_start:contig_end+1]
                front_cigar = elongate_cigar(dict_anchor_align_single["front"]["cigar"])
                tail_cigar = elongate_cigar(dict_anchor_align_single["tail"]["cigar"])[::-1]
            else:
                contig_start = dict_anchor_align_single["tail"]["ref_alignstart"]
                contig_end = dict_anchor_align_single["front"]["ref_alignend"]

                anchored_sequence = get_reverse_complement_sequence(dict_original_fasta[contigname][contig_start:contig_end+1])
                front_cigar = elongate_cigar(dict_anchor_align_single["front"]["cigar"])[::-1]
                tail_cigar = elongate_cigar(dict_anchor_align_single["tail"]["cigar"])

            dict_sequence_info[f"{regname}__{name}__{contigname}"] = {
                "sequence":anchored_sequence,
                "start":contig_start,
                "end":contig_end,
                "front_anchorlen":dict_anchor_align_single["front"]["front_anchor_alignlen"],
                "tail_anchorlen":dict_anchor_align_single["tail"]["tail_anchor_alignlen"],
                "front_anchoridy":dict_anchor_align_single["front"]["front_identity"],
                "tail_anchoridy":dict_anchor_align_single["tail"]["tail_identity"],
                "front_cigar":front_cigar,
                "tail_cigar":tail_cigar
            }
    return dict_sequence_info

def get_reverse_complement_sequence(sequence):
    dict_reverse = {
        'A' : 'T',
        'T' : 'A',
        'G' : 'C',
        'C' : 'G',
        'N' : 'N'
    }

    reversed_sequence = sequence[::-1]
    reverse_complement_sequence = list(map(dict_reverse.__getitem__, reversed_sequence))
    return ''.join(reverse_complement_sequence)

def save_anchored_sequence_with_info(dict_anchored_sequence, path_output_seq_with_anchor, path_output_seq_without_anchor, path_output_seq_info):
    dict_sequence_with_anchor = dict()
    dict_sequence_without_anchor = dict()

    file_write_info = open(path_output_seq_info, 'w')
    file_write_info.write(get_context_to_write_from_list(["region", "contig", "front_anchor_len", "tail_anchor_len", "replace_len", "front_anchor_idy", "tail_anchor_idy"]))
    
    for sequencename in dict_anchored_sequence.keys():
        ind_sequence_after_front_anchor = get_index_of_reference_matched_with_query(dict_anchored_sequence[sequencename]["front_cigar"], dict_anchored_sequence[sequencename]["front_anchorlen"])
        ind_sequence_after_tail_anchor = get_index_of_reference_matched_with_query(dict_anchored_sequence[sequencename]["tail_cigar"], dict_anchored_sequence[sequencename]["tail_anchorlen"])

        sequence_with_anchor = dict_anchored_sequence[sequencename]["sequence"]
        sequence_without_anchor = sequence_with_anchor[ind_sequence_after_front_anchor:-ind_sequence_after_tail_anchor]

        dict_sequence_with_anchor[sequencename] = sequence_with_anchor
        dict_sequence_without_anchor[sequencename] = sequence_without_anchor

        file_write_info.write(get_context_to_write_from_list([sequencename.split("__")[0], sequencename.split("__")[-1], ind_sequence_after_front_anchor, ind_sequence_after_tail_anchor, len(sequence_without_anchor), dict_anchored_sequence[sequencename]["front_anchoridy"], dict_anchored_sequence[sequencename]["tail_anchoridy"]]))
    file_write_info.close()

    save_dict_fasta_into_file(dict_sequence_with_anchor, path_output_seq_with_anchor)
    save_dict_fasta_into_file(dict_sequence_without_anchor, path_output_seq_without_anchor)

def get_context_to_write_from_list(list_values, sep = '\t', end = '\n'):
    return sep.join(list(map(str, list_values))) + end

def get_index_of_reference_matched_with_query(list_cigar_char, len_on_query):
    ind_on_ref = 0
    cnt_query_len = 0

    for cigar_char in list_cigar_char:
        if cigar_char == '=' or cigar_char == 'X':
            ind_on_ref += 1
            cnt_query_len += 1
        elif cigar_char == 'D':
            ind_on_ref += 1
        elif cigar_char == 'I':
            cnt_query_len += 1
        else:
            raise Exception(cigar_char)
        
        if cnt_query_len == len_on_query:
            break
    
    return ind_on_ref

def elongate_cigar(cigar):
    list_cigar_char = list()
    cigar_re = re.compile("([0-9]+)([A-Z,=]+)")
    for find_result in cigar_re.finditer(cigar):
        num_elong, cchar = find_result.groups()
        list_cigar_char.extend([cchar] * int(num_elong))
    return list_cigar_char


#%%
if __name__ == "__main__":
    curr_time = int(time.time())
    
    dict_path_assembly_format = {
        # Hifiasm Binned (HiFi)
        "Hifiasm_Binned_HiFi":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20/Hap{hap}/KOREF_S1.Binned_Assembly.hifi_L1000_QV20.triobin_{triobin}.haploid.l0.asm.bp.p_ctg.fasta",

        #Hifiasm Binned (Alt) (HiFi)
        "Hifiasm_Binned_ALT_HiFi":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20/Hap{hap}/KOREF_S1.Binned_Assembly.hifi_L1000_QV20.triobin_{triobin}.haploid.l0.asm.bp.a_ctg.fasta",

        # Hifiasm Binned with0 (HiFi)
        "Hifiasm_Binned_HiFi_with0":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20.with_0/Hap{hap}/KOREF_S1.Binned_Assembly.hifi_L1000_QV20.triobin_{triobin}_0.haploid.l0.asm.bp.p_ctg.fasta",
        
        # Hifiasm Binned with0 (Alt) (HiFi)
        "Hifiasm_Binned_ALT_HiFi_with0":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20.with_0/Hap{hap}/KOREF_S1.Binned_Assembly.hifi_L1000_QV20.triobin_{triobin}_0.haploid.l0.asm.bp.a_ctg.fasta",

        # Flye Binned (ONT)
        "Flye_Binned_ONT":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/ONT_UL.Filter_L50kb_QV20_Q30Perc80/Hap{hap}/assembly.fasta",
        
        # # Flye Binned with 0 (ONT)
        "Flye_Binned_ONT.with0":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/ONT_UL.Filter_L50kb_QV20_Q30Perc80.with_0/Hap{hap}/assembly.fasta",

        # Flye Binned with0 (HiFi)
        "Flye_Binned_HiFi_with0":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20.with_0/Hap{hap}/assembly.fasta",

        # Flye Binned (HiFi)
        "Flye_Binned_HiFi":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20/Hap{hap}/assembly.fasta"
    }

    dict_hap_to_parental = {
        1:'p',
        2:'m'
    }
    
    n_round = 5
    afterround = 2
    
    for hap in [1, 2]:
        for name, path_format in dict_path_assembly_format.items():
            triobin = dict_hap_to_parental[hap]

            path_original_fasta = path_format.replace("{hap}", str(hap)).replace("{triobin}", triobin)

            path_alignment_on_anchor = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Manual_Polishing/Round{n_round}_afterR{afterround}/Anchored_Sequences/Filtered_Alignments/Hap{hap}/Filter_Anchor_Alignment.on.{name}.KOREF_S1.Trio_Specific.UpperCase.Polish_R{n_round-1}.Man.hap{hap}.extra_range250.idy95.a2h95.align1k.hang1k.tsv"

            path_output_prefix = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Manual_Polishing/Round{n_round}_afterR{afterround}/Anchored_Sequences/Filtered_Sequences/Hap{hap}/Anchored_Sequence.{name}.KOREF_S1.Trio_Specific.UpperCase.Polish_R{n_round-1}.Man.hap{hap}"

            run(
                path_alignment_on_anchor,
                path_original_fasta,
                path_output_prefix,
                name
            )

   
# %%
