#%%
import os,sys
from copy import deepcopy

import pandas as pd
import numpy as np


def run(
    path_scaff_name,
    path_scaff_info,
    list_paf_name,
    list_paf_path,
    path_save,
    min_idy,
    min_align,
    max_hang,
    min_align2hang,
    max_refhang_candidate
):
    table_scaffname = pd.read_csv(path_scaff_name, sep = '\t', names = ["scaffname", "originalname"])
    dict_original2scaffname = dict(zip(table_scaffname["originalname"], table_scaffname["scaffname"].astype(str)))
    dict_scaff_info = read_scaff_info(path_scaff_info)

    table_tot_paf_aligned = read_all_paf_tables_with_filtering(list_paf_name, list_paf_path, min_idy, min_align, max_hang, min_align2hang)
    
    table_tot_paf_aligned.to_csv(path_save, sep = '\t', index = False)
    # dict_candidate_alignments = get_candidate_alignment(dict_scaff_info, table_tot_paf_aligned, max_refhang_candidate)

    # print_whether_candidate_exists(table_scaffname, dict_candidate_alignments)
    
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

def read_all_paf_tables_with_filtering(list_paf_name, list_paf_path, min_idy, min_align, max_hang, min_align2hang):
    list_paf_aligned = list()
    for name, path_paf in zip(list_paf_name, list_paf_path):
        print(name)
        assert os.path.exists(path_paf), path_paf
        table_paf_name = read_paf(path_paf)
        if table_paf_name.shape[0] == 0: 
            print("Skip Empty DataFrame: ", name)
            continue
        table_paf_name_filt_end_aligned = filter_paf_of_end_aligned_query(table_paf_name, min_idy, min_align, max_hang, min_align2hang)

        table_paf_name_filt_end_aligned.loc[:, "qryname"] = table_paf_name_filt_end_aligned["qryname"].apply(lambda val: f"{name}::{val}")
        list_paf_aligned.append(table_paf_name_filt_end_aligned)
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

def get_candidate_alignment(dict_scaff_info, table_tot_paf_aligned, max_refhang_candidate):
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
                dict_candidate_alignments[scaffname][contigind]["candidate"]["front"]["alignment"] = filter_candidate_alignment_with_direction(table_tot_paf_aligned, contigname, True, max_refhang_candidate)
            if check_end_candidate:
                dict_candidate_alignments[scaffname][contigind]["candidate"]["end"]["alignment"] = filter_candidate_alignment_with_direction(table_tot_paf_aligned, contigname, False, max_refhang_candidate)
    return dict_candidate_alignments
    
def filter_candidate_alignment_with_direction(table_paf, contigname, is_front, max_refhang_candidate):
    table_paf_of_contig = table_paf[table_paf["refname"].astype(str) == contigname]

    if is_front:
        table_paf_of_contig_direction = table_paf_of_contig[table_paf_of_contig["ref_alignstart"] < (table_paf_of_contig["reflen"]/2)]
        table_paf_of_contig_direction = table_paf_of_contig_direction[table_paf_of_contig_direction["ref_alignstart"] < max_refhang_candidate]
    else:
        table_paf_of_contig_direction = table_paf_of_contig[table_paf_of_contig["ref_alignend"] > (table_paf_of_contig["reflen"]/2)]
        table_paf_of_contig_direction = table_paf_of_contig_direction[(table_paf_of_contig_direction["reflen"] - table_paf_of_contig_direction["ref_alignend"]) < max_refhang_candidate]
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

#%%
if __name__ == "__main__":
    hap = '1'

    path_scaff_name = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon/Round13.TGS_Gapcloser_Modified_V2_re.Unpl_Scf.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.Manual_Scaff.HiFi_trio.Hifiasm_Loc.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.ONT_trio.ONT_Amb.Manual_Scaff/Prealignment_for_Candidate_Filtering/Hap{hap}/CuratedDraftASM.KOREF_S1.Trio_Specific.hap{hap}.GapFilling.Unplaced_Scf.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.Manual_Scaff.HiFi_trio.Hifiasm_Loc.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.ONT_trio.ONT_Amb.Manual_Scaff.name_map"
    path_scaff_info = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon/Round13.TGS_Gapcloser_Modified_V2_re.Unpl_Scf.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.Manual_Scaff.HiFi_trio.Hifiasm_Loc.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.ONT_trio.ONT_Amb.Manual_Scaff/Prealignment_for_Candidate_Filtering/Hap{hap}/CuratedDraftASM.KOREF_S1.Trio_Specific.hap{hap}.GapFilling.Unplaced_Scf.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.Manual_Scaff.HiFi_trio.Hifiasm_Loc.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.ONT_trio.ONT_Amb.Manual_Scaff.orignial_scaff_infos"
    list_paf_name = [
        "Hifiasm_Binned_pri",
        "Hifiasm_Binned_alt",
        "Flye_HiFi_Binned",
        "Flye_ONT_Binned",
        "Unplaced_Scaffold",
        "Flye_HiFi_Binned_with0",
        "Flye_ONT_Binned_with0",
        "Flye_HiFi_Binned_with0a",
        "Flye_ONT_Binned_with0a",
        "Flye_ONT_Binned_short"
    ]
    list_paf_path = list(map(lambda name:
                             f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon/Round13.TGS_Gapcloser_Modified_V2_re.Unpl_Scf.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.Manual_Scaff.HiFi_trio.Hifiasm_Loc.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.ONT_trio.ONT_Amb.Manual_Scaff/Prealignment_for_Candidate_Filtering/Hap{hap}/Align.Option_Adjust/Align.REF_KOREFDraftASM_Hap{hap}_Split.QRY_{name}.asm20Pk15w5e0U1r50020knolongjoins100z800100qoccfrac0f00001.filt_match100_idy0.5.paf",
                             list_paf_name))
    # list_paf_path = list(map(lambda name:
    #                          f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon/Round13.TGS_Gapcloser_Modified_V2_re.Unpl_Scf.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.Manual_Scaff.HiFi_trio.Hifiasm_Loc.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.ONT_trio.ONT_Amb.Manual_Scaff/Prealignment_for_Candidate_Filtering/Hap{hap}/Align/Align.REF_KOREFDraftASM_Hap{hap}_Split.QRY_{name}.asm5P.paf",
    #                          list_paf_name))
    for path in list_paf_path:
        assert os.path.exists(path), path

    path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon/Round13.TGS_Gapcloser_Modified_V2_re.Unpl_Scf.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.Manual_Scaff.HiFi_trio.Hifiasm_Loc.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.ONT_trio.ONT_Amb.Manual_Scaff/Prealignment_for_Candidate_Filtering/Hap{hap}/Align.Option_Adjust/Align.REF_KOREFDraftASM_Hap{hap}_Split.QRY_Any.asm20Pk15w5e0U1r50020knolongjoins100z800100qoccfrac0f00001.filtered.idy_0.7.align_5000.align2hang_0.5.paf"

    min_idy = 0.8
    min_align = 500
    max_hang = 2000
    min_align2hang = 0.8
    max_refhang_candidate = 10000000

    run(
        path_scaff_name,
        path_scaff_info,
        list_paf_name,
        list_paf_path,
        path_save,
        min_idy,
        min_align,
        max_hang,
        min_align2hang,
        max_refhang_candidate
    )
# %%
