#%%
import json, os, sys
from pathlib import Path

import pandas as pd
import numpy as np


def run(path_ragtag_agp_all, path_ragtag_paf_all, path_ragtag_agp_teloonly, path_telomere_end_info, list_chrnames_assign, manual_placed, list_scaffolds_ignore, path_save_prefix):
    os.makedirs(os.path.dirname(path_save_prefix), exist_ok=True)

    table_agp_all = pd.read_csv(path_ragtag_agp_all, sep = '\t', skiprows = 2, names = ["ref", "start", "end", "order", "type", "qry", "qry_start", "qry_end", "strand"])
    table_agp_teloonly = pd.read_csv(path_ragtag_agp_teloonly, sep = '\t', skiprows = 2, names = ["ref", "start", "end", "order", "type", "qry", "qry_start", "qry_end", "strand"])

    table_paf_all = pd.read_csv(path_ragtag_paf_all, sep = '\t', names = ["qry", "qrylen", "qrystart", "qryend", "strand", "ref", "reflen", "refstart", "refend", "matchlen", "alignlen", "quality"], usecols = list(range(12)))
    
    table_telo_info = pd.read_csv(path_telomere_end_info, sep = '\t')
    dict_scaffold_to_teloend = dict(zip(table_telo_info["Scaffold"], table_telo_info["Telomere_End"]))
    
    path_save_telo_end_scaffolds = f"{path_save_prefix}.telo_end_assigned.json"
    dict_telo_end_scaffold = assign_telomere_ends(table_agp_teloonly, dict_scaffold_to_teloend, list_chrnames_assign, manual_placed, path_save_telo_end_scaffolds)

    path_save_reassigned_order = f"{path_save_prefix}.reordered_scaffolds.json"
    dict_reassigned_order = reassign_ragtag_order_with_telomere_assignment(table_agp_all, table_paf_all, list_chrnames_assign, list_scaffolds_ignore, dict_telo_end_scaffold, path_save_reassigned_order)

def assign_telomere_ends(table_agp_teloonly, dict_scaffold_to_teloend, list_chrnames_assign, manual_placed, path_save):
    table_agp_for_assign = table_agp_teloonly[table_agp_teloonly["ref"].apply(lambda val: val.split('_')[0] in list_chrnames_assign)]

    table_agp_for_assign_nongap = table_agp_for_assign[~table_agp_for_assign["type"].isin(['U', 'N'])]

    dict_telo_assigned = dict()

    for refname in set(list(table_agp_for_assign_nongap["ref"].unique())+list(manual_placed.keys())):
        refname_simplify = refname.split('_')[0]
        table_agp_for_assign_singlechr = table_agp_for_assign_nongap[table_agp_for_assign_nongap["ref"] == refname]

        dict_telo_tmp = {
            "front":list(),
            "end":list()
        }

        for _, row in table_agp_for_assign_singlechr.iterrows():
            scfname = row["qry"]
            telo_direction = get_telo_direction_from_scfname(scfname, dict_scaffold_to_teloend, row["strand"])
            dict_telo_tmp[telo_direction].append((scfname, row["order"], row["strand"]))
        
        assert len(dict_telo_tmp["front"]) < 2, dict_telo_tmp
        assert len(dict_telo_tmp["end"]) < 2, dict_telo_tmp
        
        if len(dict_telo_tmp["front"]) == 1 and len(dict_telo_tmp["end"]) == 1:
            assert dict_telo_tmp["front"][0][1] < dict_telo_tmp["end"][0][1]
        
        if manual_placed.get(refname_simplify) != None:
            for manual_direction, manual_scfname in manual_placed[refname_simplify]:
                assert len(dict_telo_tmp.get(manual_direction)) == 0, dict_telo_tmp
                manual_scfname_telo_direction = get_telo_direction_from_scfname(manual_scfname, dict_scaffold_to_teloend, '+')
                if manual_direction != manual_scfname_telo_direction:
                    manual_scfname_direction = '-'
                    print(manual_scfname, manual_direction, manual_scfname_telo_direction)
                else:
                    manual_scfname_direction = '+'

                dict_telo_tmp[manual_direction] = [(manual_scfname, -1, manual_scfname_direction)]
        
        dict_telo_assigned[refname_simplify] = {
            "front" : (dict_telo_tmp["front"][0][0], dict_telo_tmp["front"][0][2]) if len(dict_telo_tmp["front"]) == 1 else (None, None),
            "end" : (dict_telo_tmp["end"][0][0], dict_telo_tmp["end"][0][2]) if len(dict_telo_tmp["end"]) == 1 else (None, None)
        }
        print(dict_telo_tmp)
    
    with open(path_save, 'w') as fw:
        json.dump(dict_telo_assigned, fw, indent = 4)
    return dict_telo_assigned

def get_telo_direction_from_scfname(scfname, dict_telo_end, ori):
    telo_direction = dict_telo_end.get(scfname)
    if telo_direction == None:
        if "ntLink" in scfname:
            ntlink_telo_direction = scfname.split('.')[1]
            assert ntlink_telo_direction in {"Left", "Right"}, scfname
            telo_direction = "front" if ntlink_telo_direction == "Left" else "end"
        else:
            raise Exception(scfname)
    if ori == '-':
        telo_direction = {"front":"end", "end":"front"}[telo_direction]
    return telo_direction

def reassign_ragtag_order_with_telomere_assignment(table_agp, table_paf, list_chrnames_assign, list_scaffolds_ignore, dict_telo_end_scaffold, path_save_order):
    table_agp_nongap = table_agp[~table_agp["type"].isin(['U', 'N'])]
    table_agp_nongap["ref"] = table_agp_nongap["ref"].apply(lambda val: val.split('_')[0])
    dict_reassigned_order = dict()

    table_agp_nongap = table_agp_nongap[~table_agp_nongap["qry"].isin(list_scaffolds_ignore)]
    
    # table_paf = table_paf[table_paf["alignlen"] > 1000]

    for chrname in list_chrnames_assign:
        dict_reassigned_order[chrname] = list()
        # dict_reassigned_order[chrname].append(dict_telo_end_scaffold[chrname]["front"])        

        table_agp_chrname = table_agp_nongap[table_agp_nongap["ref"] == chrname]
        
        list_scaffold_in_order = list()

        is_telo_started = False
        is_telo_ended = False

        table_agp_chrname_sorted = table_agp_chrname.sort_values(by = ["order"])

        if dict_telo_end_scaffold[chrname]["front"][0] == None:
            is_telo_started = True
            list_scaffold_in_order.append((None, None, 0, 0))
        elif dict_telo_end_scaffold[chrname]["front"][0] not in table_agp_chrname_sorted["qry"].to_list():
            is_telo_started = True
            list_scaffold_in_order.append((dict_telo_end_scaffold[chrname]["front"][0], dict_telo_end_scaffold[chrname]["front"][1], 0, 0))

        for _, row in table_agp_chrname_sorted.iterrows():
            scfname = row["qry"]
            scf_ori = row["strand"]

            if scfname == dict_telo_end_scaffold[chrname]["front"][0]:
                assert scf_ori == dict_telo_end_scaffold[chrname]["front"][1]
                is_telo_started = True
            if not is_telo_started:
                continue

            # pos_align_start, pos_align_end = get_aligned_region_with_percentile(table_paf, chrname, scfname)
            table_paf_comb = table_paf[np.logical_and(table_paf["ref"] == chrname, table_paf["qry"] == scfname)]
            pos_align_start = table_paf_comb["refstart"].min()
            pos_align_end = table_paf_comb["refend"].max()
            
            list_scaffold_in_order.append((scfname, scf_ori, pos_align_start, pos_align_end))

            if scfname == dict_telo_end_scaffold[chrname]["end"][0]:
                assert scf_ori == dict_telo_end_scaffold[chrname]["end"][1]
                is_telo_ended = True
            if is_telo_ended:
                break
        
        if dict_telo_end_scaffold[chrname]["end"][0] == None:
            is_telo_ended = True
            list_scaffold_in_order.append((None, None, np.inf, np.inf))
        elif dict_telo_end_scaffold[chrname]["end"][0] not in table_agp_chrname_sorted["qry"].to_list():
            is_telo_ended = True
            list_scaffold_in_order.append((dict_telo_end_scaffold[chrname]["end"][0], dict_telo_end_scaffold[chrname]["end"][1], np.inf, np.inf))
        
        bef_start, bef_end = None, None
        for _, _, start, end in list_scaffold_in_order:
            if bef_start != None:
                # assert start > bef_end, list_scaffold_in_order
                assert start > bef_start and end > bef_end, list_scaffold_in_order
            bef_start = start
            bef_end = end
        assert check_is_sorted(list(map(lambda val: val[2], list_scaffold_in_order))), list_scaffold_in_order
        assert check_is_sorted(list(map(lambda val: val[3], list_scaffold_in_order))), list_scaffold_in_order
        
        for scfname, ori, _, _ in list_scaffold_in_order:
            dict_reassigned_order[chrname].append([scfname, ori])
        
        # dict_reassigned_order[chrname].append(dict_telo_end_scaffold[chrname]["end"])
    with open(path_save_order, 'w') as fw:
        json.dump(dict_reassigned_order, fw, indent = 4)

    return dict_reassigned_order

def get_aligned_region_with_percentile(table_paf, refname, scfname):
    from copy import deepcopy
    table_paf_of_this_comb = table_paf[np.logical_and(table_paf["qry"] == scfname, table_paf["ref"] == refname)]

    list_mapped_ranges = list()
    for _, row in table_paf_of_this_comb.iterrows():
        aligned_ref_start = row["refstart"]
        aligned_ref_end = row["refend"]

        is_range_merged = False
        for ind, ranges in enumerate(deepcopy(list_mapped_ranges)):
            rs, re = ranges

            if (aligned_ref_start >= rs and aligned_ref_start <= re) or (aligned_ref_end >= rs and aligned_ref_end <= re):
                is_range_merged = True
                list_mapped_ranges[ind] = (min(rs, aligned_ref_start), max(re, aligned_ref_end))
                break
            else:
                pass
        if not is_range_merged:
            list_mapped_ranges.append((aligned_ref_start, aligned_ref_end))
    
    uniq_positions = list()
    for ranges in list_mapped_ranges:
        uniq_positions.extend(range(ranges[0], ranges[1]))
    uniq_positions = list(set(uniq_positions))

    pos_q_front = np.quantile(uniq_positions, 0.025)
    pos_q_end = np.quantile(uniq_positions, 0.975)

    return (pos_q_front, pos_q_end)

def check_is_sorted(values):
    if len(values) == 0:
        return True
    
    bef_val = values[0]
    for val in values[1:]:
        if val < bef_val:
            return False
        bef_val = val
    return True


#%%
if __name__ == "__main__":
    path_ragtag_agp_all = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling.RagTag.default_option.CHM13/DraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.no_telo.telo_single.hap1.RagTag.default_option.minimap2_latest.CHM13/ragtag.scaffold.agp"
    path_ragtag_paf_all = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling.RagTag.default_option.CHM13/DraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.no_telo.telo_single.hap1.RagTag.default_option.minimap2_latest.CHM13/ragtag.scaffold.asm.paf"
    path_ragtag_agp_teloonly = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling.RagTag.default_option.CHM13/DraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.telo_single.hap1.RagTag.default_option.minimap2_latest.CHM13/ragtag.scaffold.agp"
    path_telomere_end_info = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/Telomere_Detection/Fasta_Split_by_Telomere/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap1.p_ctg.telo_single.telo_info.tsv"

    list_chrnames_assign = [
        "chr4", "chr13", "chr15", "chr16", "chr17", "chr18", "chr22", "chrY"
    ]
    manual_placed = {
        
    }
    list_scaffolds_ignore = []

    path_save_prefix = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling.RagTag.default_option.CHM13/CuratedDraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.non_t2t_scaffold.hap1.RagTag.default_option.minimap2_latest.CHM13/CuratedDraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.hap1.RagTag_with_Telomere_Assign.minimap2_latest.CHM13"

    run(path_ragtag_agp_all, path_ragtag_paf_all, path_ragtag_agp_teloonly, path_telomere_end_info, list_chrnames_assign, manual_placed, list_scaffolds_ignore, path_save_prefix)

#%%
if __name__ == "__main__":
    path_ragtag_agp_all = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling.RagTag.default_option.CHM13/DraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.no_telo.telo_single.hap2.RagTag.default_option.minimap2_latest.CHM13/ragtag.scaffold.agp"
    path_ragtag_paf_all = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling.RagTag.default_option.CHM13/DraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.no_telo.telo_single.hap2.RagTag.default_option.minimap2_latest.CHM13/ragtag.scaffold.asm.paf"
    path_ragtag_agp_teloonly = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling.RagTag.default_option.CHM13/DraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.telo_single.hap2.RagTag.default_option.minimap2_latest.CHM13/ragtag.scaffold.agp"
    path_telomere_end_info = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/Telomere_Detection/Fasta_Split_by_Telomere/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap2.p_ctg.telo_single.telo_info.tsv"

    
    list_chrnames_assign = [
        "chr4", "chr14", "chr15", "chr17", "chr18", "chr19", "chr22"
    ]
    manual_placed = {
        "chr15":[("front", "h2tg000033l")]
    }
    list_scaffolds_ignore = ["h2tg000045l"]

    path_save_prefix = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling.RagTag.default_option.CHM13/CuratedDraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.non_t2t_scaffold.hap2.RagTag.default_option.minimap2_latest.CHM13/CuratedDraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.hap2.RagTag_with_Telomere_Assign.minimap2_latest.CHM13"

    run(path_ragtag_agp_all, path_ragtag_paf_all, path_ragtag_agp_teloonly, path_telomere_end_info, list_chrnames_assign, manual_placed, list_scaffolds_ignore, path_save_prefix)
# %%
