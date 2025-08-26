#%%
import os,sys
from copy import deepcopy
from joblib import Parallel, delayed

import pandas as pd
import numpy as np


def run(
    path_scaff_name,
    path_scaff_info,
    list_paf_name,
    list_paf_path,
    min_idy,
    min_align,
    max_hang,
    min_align2hang,
    max_refhang_candidate,
    min_elongate,
    path_output_prefix,
    n_jobs = 10
):
    table_scaffname = pd.read_csv(path_scaff_name, sep = '\t', names = ["scaffname", "originalname"])
    dict_original2scaffname = dict(zip(table_scaffname["originalname"], table_scaffname["scaffname"].astype(str)))
    dict_scaff_info = read_scaff_info(path_scaff_info)

    table_tot_paf_aligned = read_all_paf_tables_with_filtering(list_paf_name, list_paf_path, min_idy, min_align, max_hang, min_align2hang, n_jobs = n_jobs)
    
    dict_candidate_alignments = get_candidate_alignment(dict_scaff_info, table_tot_paf_aligned, max_refhang_candidate, min_elongate)

    print_whether_candidate_exists(table_scaffname, dict_candidate_alignments)
    
    list_scaffnames_with_gap = list(filter(lambda name: len(dict_scaff_info[name]) > 1, dict_scaff_info.keys()))
    list_original_chrnames_with_gap = list(filter(lambda chrname: dict_original2scaffname[chrname] in list_scaffnames_with_gap, dict_original2scaffname.keys()))
    
    for chrname in list_original_chrnames_with_gap:
        save_candidate_anchors_as_table(dict_candidate_alignments, chrname, dict_original2scaffname, path_output_prefix)
    
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

def read_all_paf_tables_with_filtering(list_paf_name, list_paf_path, min_idy, min_align, max_hang, min_align2hang, n_jobs = 1):
    list_paf_aligned = list()
    
    def read_paf_with_filtering(path_paf, name_paf):
        assert os.path.exists(path_paf), path_paf
        table_paf_name = read_paf(path_paf)
        table_paf_name_filt_end_aligned = filter_paf_of_end_aligned_query(table_paf_name, min_idy, min_align, max_hang, min_align2hang)

        table_paf_name_filt_end_aligned.loc[:, "qryname"] = table_paf_name_filt_end_aligned["qryname"].apply(lambda val: f"{name_paf}::{val}")
        return table_paf_name_filt_end_aligned
    
    with Parallel(n_jobs = n_jobs) as parallel:
        list_paf_aligned = parallel(delayed(read_paf_with_filtering)(
            path_paf, name
        )for (name, path_paf) in zip(list_paf_name, list_paf_path))
        
    table_tot_paf_aligned = pd.concat(list_paf_aligned)
    return table_tot_paf_aligned

def filter_paf_of_end_aligned_query(table_paf, min_idy, min_align, max_hang, min_align2hang):
    table_paf["identity"] = table_paf["match"] / table_paf["align"]

    table_paf_filt = table_paf[table_paf["identity"] > min_idy]
    table_paf_filt = table_paf_filt[table_paf_filt["align"] > min_align]

    table_paf_filt["wronghang"] = table_paf_filt.apply(get_wrong_hang_length, axis = 1)
    table_paf_filt["aligned2hang"] = table_paf_filt["align"] / (table_paf_filt["align"]+table_paf_filt["wronghang"])

    table_paf_filt = table_paf_filt[table_paf_filt["wronghang"] < max_hang]
    table_paf_filt = table_paf_filt[table_paf_filt["aligned2hang"] > min_align2hang]

    return table_paf_filt
    # return table_paf

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

def get_candidate_alignment(dict_scaff_info, table_tot_paf_aligned, max_refhang_candidate, min_elongate):
    dict_candidate_alignments = deepcopy(dict_scaff_info)

    for scaffname in dict_scaff_info.keys():
        list_contigind = list(dict_scaff_info[scaffname].keys())
        for contigind in list_contigind:
            check_front_candidate = True
            check_end_candidate = True

            if contigind == list_contigind[0]:
                check_front_candidate = False
            if contigind == list_contigind[-1]:
                check_end_candidate = False
            
            contigname = dict_scaff_info[scaffname][contigind]["name"]

            dict_candidate_alignments[scaffname][contigind]["candidate"] = {
                "front" : {
                    "required" : check_front_candidate,
                    "alignment" : dict()
                },
                "end" : {
                    "required" : check_end_candidate,
                    "alignment" : dict()
                }
            }
            if check_front_candidate:
                dict_candidate_alignments[scaffname][contigind]["candidate"]["front"]["alignment"] = filter_candidate_alignment_with_direction(table_tot_paf_aligned, contigname, True, max_refhang_candidate, min_elongate)
            if check_end_candidate:
                dict_candidate_alignments[scaffname][contigind]["candidate"]["end"]["alignment"] = filter_candidate_alignment_with_direction(table_tot_paf_aligned, contigname, False, max_refhang_candidate, min_elongate)
    return dict_candidate_alignments
    
def filter_candidate_alignment_with_direction(table_paf, contigname, is_front, max_refhang_candidate, min_elongate):
    table_paf_of_contig = table_paf[table_paf["refname"].astype(str) == contigname]

    if is_front:
        table_paf_of_contig_direction = table_paf_of_contig[table_paf_of_contig["ref_alignstart"] < (table_paf_of_contig["reflen"]/2)]
        table_paf_of_contig_direction = table_paf_of_contig_direction[table_paf_of_contig_direction["ref_alignstart"] < max_refhang_candidate]

        table_paf_of_contig_direction["elongated_length"] = table_paf_of_contig_direction.apply(lambda row: row["qry_alignstart"] if row["strand"] == '+' else row["qrylen"]-row["qry_alignend"], axis = 1)
        table_paf_of_contig_direction = table_paf_of_contig_direction[table_paf_of_contig_direction["elongated_length"] > min_elongate]

    else:
        table_paf_of_contig_direction = table_paf_of_contig[table_paf_of_contig["ref_alignend"] > (table_paf_of_contig["reflen"]/2)]
        table_paf_of_contig_direction = table_paf_of_contig_direction[(table_paf_of_contig_direction["reflen"] - table_paf_of_contig_direction["ref_alignend"]) < max_refhang_candidate]

        table_paf_of_contig_direction["elongated_length"] = table_paf_of_contig_direction.apply(lambda row: row["qry_alignstart"] if row["strand"] == '-' else row["qrylen"]-row["qry_alignend"], axis = 1)
        table_paf_of_contig_direction = table_paf_of_contig_direction[table_paf_of_contig_direction["elongated_length"] > min_elongate]

    return table_paf_of_contig_direction

def print_whether_candidate_exists(table_scaffname, dict_candidate_alignments):
    table_scaffname.loc[:, "scaffname"] = table_scaffname["scaffname"].astype(str)
    dict_scaffname_to_originalname = dict(zip(table_scaffname["scaffname"], table_scaffname["originalname"]))
    for scaffname in dict_scaffname_to_originalname.keys():
        if len(dict_candidate_alignments[scaffname]) == 1:
            # print(f"{dict_scaffname_to_originalname[scaffname]} Finished", end = "\n\n")
            pass
        else:
            print(f"{dict_scaffname_to_originalname[scaffname]}")
            for contigind in dict_candidate_alignments[scaffname].keys():
                contigname = dict_candidate_alignments[scaffname][contigind]["name"]
                print(contigname)
                if dict_candidate_alignments[scaffname][contigind]["candidate"]["front"]["required"]:
                    n_tot_aligns = dict_candidate_alignments[scaffname][contigind]["candidate"]["front"]["alignment"].shape[0]
                    print(f"\tFront\t{n_tot_aligns}")
                if dict_candidate_alignments[scaffname][contigind]["candidate"]["end"]["required"]:
                    n_tot_aligns = dict_candidate_alignments[scaffname][contigind]["candidate"]["end"]["alignment"].shape[0]
                    print(f"\tEnd\t{n_tot_aligns}")
            print()

def save_candidate_anchors_as_table(dict_candidate_alignments, chrname, dict_original2scaffname, path_output_prefix):
    for contigind, dict_aligned in dict_candidate_alignments[dict_original2scaffname[chrname]].items():
        if dict_aligned["candidate"]["front"]["required"]:
            table_candidate_anchors = dict_aligned["candidate"]["front"]["alignment"]
            path_save = f"{path_output_prefix}.chrname_{chrname}.contig_index_{contigind}.Front.tsv"
            table_candidate_anchors.to_csv(path_save, sep = '\t', index = False)

        if dict_aligned["candidate"]["end"]["required"]:
            table_candidate_anchors = dict_aligned["candidate"]["end"]["alignment"]
            path_save = f"{path_output_prefix}.chrname_{chrname}.contig_index_{contigind}.End.tsv"
            table_candidate_anchors.to_csv(path_save, sep = '\t', index = False)


#%%
if __name__ == "__main__":
    hap = '1'

    path_scaff_name = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon_wo_RawONT/Round8.Manual_Scaff/Hap{hap}/Prealignment/tgsseqsplit/CuratedDraftASM.KOREF_S1.Trio_Specific.hap{hap}.Round8.GapFilling.Manual_Scaff.LocalASM.Unplaced_Scf.Hifiasm_Hifi_Bin.Hifiasm_Hifi_Bin0.Flye_Hifi_Bin.Flye_Hifi_Bin0.Manual_Scaff.name_map"
    path_scaff_info = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon_wo_RawONT/Round8.Manual_Scaff/Hap{hap}/Prealignment/tgsseqsplit/CuratedDraftASM.KOREF_S1.Trio_Specific.hap{hap}.Round8.GapFilling.Manual_Scaff.LocalASM.Unplaced_Scf.Hifiasm_Hifi_Bin.Hifiasm_Hifi_Bin0.Flye_Hifi_Bin.Flye_Hifi_Bin0.Manual_Scaff.orignial_scaff_infos"
    list_paf_name = [
        "Hifiasm_Binned_HiFi",
        "Hifiasm_Binned_ALT_HiFi",
        "Hifiasm_Binned_HiFi_with0",
        "Hifiasm_Binned_ALT_HiFi_with0",
        "Flye_Binned_HiFi_with0",
        "Flye_Binned_HiFi",
        "Unplaced_Scaffold"
    ]
    list_paf_path = list(map(lambda name:
                             f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon_wo_RawONT/Round8.Manual_Scaff/Hap{hap}/Prealignment/Align.Original_Contig/Align.REF_KOREFDraftASM_Hap{hap}_Split.QRY_{name}.asm5Pc.filt_match100_idy0.2.paf",
                             list_paf_name))


    min_idy = 0.8
    min_align = 500
    max_hang = 10000
    min_align2hang = 0.9
    max_refhang_candidate = 10000
    min_elongate = 1000

    path_output_prefix = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon_wo_RawONT/Round8.Manual_Scaff/Hap{hap}/Prealignment/Anchored_Alignments/Anchored_Alignments.REF_KOREFDraftASM_Hap{hap}_Split"

    os.makedirs(os.path.dirname(path_output_prefix), exist_ok = True)
    # run(
    #     path_scaff_info,
    #     list_paf_name,
    #     list_paf_path,
    #     min_idy,
    #     min_align,
    #     max_hang
    # )
# %%
