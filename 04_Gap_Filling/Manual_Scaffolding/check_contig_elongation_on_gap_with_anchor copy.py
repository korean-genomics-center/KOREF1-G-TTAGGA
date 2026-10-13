#%%
from pathlib import Path
from joblib import Parallel, delayed
from collections import Counter
import sys,os

import pandas as pd
import numpy as np

sys.path.append(str(Path(__file__).parents[0]))
from check_anchored_contig import read_all_paf_tables_with_filtering, get_candidate_alignment
from check_contig_overlap_elongation import read_all_pairwise_paf_tables_with_filtering, generate_contig_linkage_by_overlap, get_shortest_contig_paths_from_specific_contig_with_single_direction, get_statistics_of_elongated_path_on_paths_as_table, get_ranges_for_scaffold_elongation_from_path


def main(
    path_scaff_name,
    path_scaff_info,
    list_paf_name_against_ref,
    list_paf_path_against_ref,
    list_paf_name_pairwise,
    list_paf_path_pairwise,
    min_idy = 0.8,
    min_align_contig_to_contig = 500,
    min_align_contig_to_ref = 5000,
    max_hang_contig_to_contig = 10000,
    max_hang_contig_to_ref = 1000000,
    min_align2hang_contig_to_contig = 0.8,
    min_align2hang_contig_to_ref = 0.5,
    n_alignedcontig_check = 50,
    max_distance_ratio_from_last_contig = 0.3,
    n_jobs = 20
):
    table_scaffname = pd.read_csv(path_scaff_name, sep = '\t', names = ["scaffname", "originalname"])
    dict_original2scaffname = dict(zip(table_scaffname["originalname"], table_scaffname["scaffname"].astype(str)))
    dict_scaff_info = read_scaff_info(path_scaff_info)
    list_contignames_between_gap = get_contig_pairs_with_gap(dict_scaff_info)

    table_all_refaligned_paf_filt = read_all_paf_tables_with_filtering(list_paf_name_against_ref, list_paf_path_against_ref, min_idy, min_align_contig_to_ref, max_hang_contig_to_ref, min_align2hang_contig_to_ref)
    # list_contignames_aligned_multiple_against_ref = table_all_refaligned_paf_filt[table_all_refaligned_paf_filt["qryname"].duplicated()]["qryname"].unique()
    # table_all_refaligned_paf_filt_unique = table_all_refaligned_paf_filt[~table_all_refaligned_paf_filt["qryname"].isin(list_contignames_aligned_multiple_against_ref)]
    table_all_refaligned_paf_filt["refname"] = table_all_refaligned_paf_filt["refname"].astype(str)

    table_all_pairwise_paf_filt = read_all_pairwise_paf_tables_with_filtering(list_paf_name_pairwise, list_paf_path_pairwise, min_idy, min_align_contig_to_contig, max_hang_contig_to_contig, min_align2hang_contig_to_contig)
    dict_contig_pairwise_alignment = generate_contig_linkage_by_overlap(table_all_pairwise_paf_filt)

    dict_chain_stats = find_elongation_path_from_aligned_contigs(list_contignames_between_gap, table_all_refaligned_paf_filt, dict_contig_pairwise_alignment, n_alignedcontig_check, max_distance_ratio_from_last_contig, n_jobs)

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

def get_contig_pairs_with_gap(dict_scaff_info):
    list_contig_pairs = list()

    for scaffname, dict_contigs in dict_scaff_info.items():
        list_contig_order = list(dict_contigs.keys())
        if len(list_contig_order) > 1:
            for ind1, contigorder1 in enumerate(list_contig_order[:-1]):
                for contigorder2 in list_contig_order[ind1+1:]:
                    # contigorder2 = list_contig_order[ind1+1]
                    list_contig_pairs.append(
                        (dict_contigs[contigorder1]["name"], dict_contigs[contigorder2]["name"])
                    )
    return list_contig_pairs

def find_elongation_path_from_aligned_contigs(list_contignames_between_gap, table_all_refaligned_paf_filt, dict_contig_pairwise_alignment, n_alignedcontig_check, max_distance_ratio_from_last_contig, n_jobs):
    dict_elongation_path = dict()

    for contignames in list_contignames_between_gap:
        contignames_key = "::".join(contignames)
        dict_elongation_path[contignames_key] = pd.DataFrame()

        table_contigs_aligned_depart = get_aligned_contigs_against_ref_from_end(table_all_refaligned_paf_filt, contignames[0], False, n_alignedcontig_check, max_distance_ratio_from_last_contig)
        table_contigs_aligned_dest = get_aligned_contigs_against_ref_from_end(table_all_refaligned_paf_filt, contignames[1], True, n_alignedcontig_check, max_distance_ratio_from_last_contig)

        print(contignames_key, table_contigs_aligned_depart.shape[0], table_contigs_aligned_dest.shape[0])
        if table_contigs_aligned_depart.shape[0] == 0 or table_contigs_aligned_dest.shape[0] == 0:
            continue
        
        list_depart_contig_to_strand = list(table_contigs_aligned_depart.apply(lambda row: f"{row['qryname']}__strand__{row['strand']}", axis = 1).unique())
        dict_contig_dest_to_strand = dict()
        for _, row in table_contigs_aligned_dest.drop_duplicates(subset = ["qryname", "strand"]).iterrows():
            qryname_dest = row["qryname"]
            strand_dest = row["strand"]
            if dict_contig_dest_to_strand.get(qryname_dest) == None: dict_contig_dest_to_strand[qryname_dest] = list()
            dict_contig_dest_to_strand[qryname_dest].append(strand_dest)

        dict_contig_path_from_depart_to_dest = get_paths_from_depart_to_dest_contig_parallel(list_depart_contig_to_strand, dict_contig_dest_to_strand, dict_contig_pairwise_alignment, n_jobs)
        # dict_contig_path_from_depart_to_dest = get_paths_from_depart_to_dest_contig_parallel(dict_contig_depart_to_strand, dict_contig_dest_to_strand, dict_contig_pairwise_alignment, n_jobs)
        table_path_stats = convert_depart_to_dest_paths_into_stats(dict_contig_path_from_depart_to_dest, list_depart_contig_to_strand, dict_contig_pairwise_alignment)
        dict_elongation_path[contignames_key] = table_path_stats
    return dict_elongation_path
        
def get_aligned_contigs_against_ref_from_end(table_paf, ref_contigname, get_front, n_collect, distance_limit_ratio):
    table_paf_of_thisref = table_paf[table_paf["refname"] == ref_contigname]
    if table_paf_of_thisref.shape[0] == 0:
        return table_paf_of_thisref
    
    # list_qrynames_amb = get_ambiguous_strand_query_from_paf(table_paf_of_thisref)
    # table_paf_of_thisref = table_paf_of_thisref[~table_paf_of_thisref["qryname"].isin(list_qrynames_amb)]

    thisref_len = table_paf_of_thisref["reflen"].to_list()[0]
    distance_limit = thisref_len * distance_limit_ratio    

    if get_front:
        table_paf_of_thisref_limited = table_paf_of_thisref[table_paf_of_thisref["ref_alignstart"] < distance_limit]
        table_paf_of_thisref_sorted_limited = table_paf_of_thisref_limited.sort_values(by = "ref_alignstart", ascending = get_front).head(n= (n_collect if n_collect else table_paf_of_thisref_limited.shape[0]))
    else:
        table_paf_of_thisref_limited = table_paf_of_thisref[table_paf_of_thisref["ref_alignend"] > (thisref_len-distance_limit)]
        table_paf_of_thisref_sorted_limited = table_paf_of_thisref_limited.sort_values(by = "ref_alignend", ascending = get_front).head(n=(n_collect if n_collect else table_paf_of_thisref_limited.shape[0]))

    return table_paf_of_thisref_sorted_limited

def get_ambiguous_strand_query_from_paf(table_paf):
    table_paf_qrystrand_dropdup = table_paf.drop_duplicates(subset = ["qryname", "strand"])
    count_qryname = Counter(table_paf_qrystrand_dropdup["qryname"])
    list_qrynames_amb = list(filter(lambda name: count_qryname[name] > 1, count_qryname.keys()))
    return list_qrynames_amb

def get_paths_from_depart_to_dest_contig_parallel(list_depart_contig_to_strand, dict_contig_dest_to_strand, dict_contig_alignment, n_jobs):
    # list_contignames_depart = list(dict_contig_depart_to_strand.keys())
    list_contignames_dest = list(dict_contig_dest_to_strand.keys())

    with Parallel(n_jobs = n_jobs) as parallel:
        list_dict_paths_until_target = parallel(
            delayed(get_paths_until_another_aligned_contig_from_specific_contig)(
                contigname_depart__strand.split("__strand__")[0],
                list_contignames_dest,
                dict_contig_alignment,
                contigname_depart__strand.split("__strand__")[1],
                dict_contig_dest_to_strand
            ) for contigname_depart__strand in list_depart_contig_to_strand
        )
    
    dict_path_between_depart_and_dest = dict()
    for contigname_depart__strand, dict_paths_until_dest in zip(list_depart_contig_to_strand, list_dict_paths_until_target):
        if len(dict_paths_until_dest.keys()) > 0:
            dict_path_between_depart_and_dest[contigname_depart__strand] = dict_paths_until_dest
    return dict_path_between_depart_and_dest

def get_paths_until_another_aligned_contig_from_specific_contig(contigname_start, list_contignames_end_target, dict_contig_alignment, contig_start_aligned_strand, dict_contig_dest_to_strand):
    dict_path_from_start = get_shortest_contig_paths_from_specific_contig_with_single_direction(dict_contig_alignment, contigname_start, contig_start_aligned_strand == '-', max_dist = 50)

    strand_exchange = {
        '+' : '-',
        '-' : '+'
    }
    dict_paths_until_target = dict()
    set_contignames_end_target = set(list_contignames_end_target)
    
    for contig_end, list_of_list_paths in dict_path_from_start.items():
        if contig_end in set_contignames_end_target:
            list_of_list_paths_with_correct_strand = list()

            for list_path in list_of_list_paths:
                dict_trimming_info_of_path = get_ranges_for_scaffold_elongation_from_path(list_path, dict_contig_alignment, "front" if contig_start_aligned_strand == '-' else "end")
                end_contig_strand = dict_trimming_info_of_path[contig_end]["strand"]
                if contig_start_aligned_strand == '-':
                    end_contig_strand = strand_exchange[end_contig_strand]
                
                if end_contig_strand in dict_contig_dest_to_strand[contig_end]:
                    list_of_list_paths_with_correct_strand.append(list_path)
                else:
                    print(contigname_start, contig_end, end_contig_strand, dict_contig_dest_to_strand[contig_end], flush = True)

            dict_paths_until_target[contig_end] = list_of_list_paths_with_correct_strand
    return dict_paths_until_target

def convert_depart_to_dest_paths_into_stats(dict_contig_path_from_depart_to_dest, list_depart_contig_to_strand, dict_contig_pairwise_alignment):
    table_path_stats = pd.DataFrame()

    for contig_depart_strand in list_depart_contig_to_strand:
        contig_depart, strand_depart = contig_depart_strand.split("__strand__")
        for contig_dest in dict_contig_path_from_depart_to_dest[contig_depart_strand].keys():
            list_paths_between_contig = dict_contig_path_from_depart_to_dest[contig_depart_strand][contig_dest]
            if len(list_paths_between_contig) > 0:
                table_elongation_stat = get_statistics_of_elongated_path_on_paths_as_table(list_paths_between_contig, dict_contig_pairwise_alignment, "end" if strand_depart == '+' else "front")
                table_elongation_stat["Depart"] = contig_depart
                table_elongation_stat["Destination"] = contig_dest
                table_elongation_stat["Length"] = table_elongation_stat[["Path", "Direction"]].apply(
                    lambda row: get_total_path_length(row["Path"].split("->"), dict_contig_pairwise_alignment, row["Direction"]),
                    axis = 1
                )
                table_path_stats = pd.concat([table_path_stats, table_elongation_stat])
    return table_path_stats

def get_total_path_length(path, dict_contig_alignment, direction):
    dict_ranges_for_trimming = get_ranges_for_scaffold_elongation_from_path(path, dict_contig_alignment, direction)
    list_all_ranges = list(map(lambda trimming_config: trimming_config["range"], dict_ranges_for_trimming.values()))
    list_range_length = list(map(lambda ranges: ranges[1]-ranges[0]+1, list_all_ranges))

    return sum(list_range_length)


#%%
if __name__ == "__main__":
    list_paf_name = [
        "Hifiasm_Binned_pri",
        "Hifiasm_Binned_alt",
        "Flye_HiFi_Binned",
        "Flye_ONT_Binned",
        "Unplaced_Scaffold",
        "Flye_HiFi_Binned_with0",
        "Flye_ONT_Binned_with0",
        "Flye_HiFi_Binned_with0a",
        "Flye_ONT_Binned_with0a"
    ]

    hap = '1'
    triobin = 'p'
    
    path_scaff_name = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon/Round13.TGS_Gapcloser_Modified_V2_re.Unpl_Scf.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.Manual_Scaff.HiFi_trio.Hifiasm_Loc.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.ONT_trio.ONT_Amb.Manual_Scaff/Prealignment_for_Candidate_Filtering/Hap{hap}/CuratedDraftASM.KOREF_S1.Trio_Specific.hap{hap}.GapFilling.Unplaced_Scf.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.Manual_Scaff.HiFi_trio.Hifiasm_Loc.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.ONT_trio.ONT_Amb.Manual_Scaff.name_map"
    path_scaff_info = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon/Round13.TGS_Gapcloser_Modified_V2_re.Unpl_Scf.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.Manual_Scaff.HiFi_trio.Hifiasm_Loc.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.ONT_trio.ONT_Amb.Manual_Scaff/Prealignment_for_Candidate_Filtering/Hap{hap}/CuratedDraftASM.KOREF_S1.Trio_Specific.hap{hap}.GapFilling.Unplaced_Scf.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.Manual_Scaff.HiFi_trio.Hifiasm_Loc.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.ONT_trio.ONT_Amb.Manual_Scaff.orignial_scaff_infos"

    list_paf_name_against_ref = list_paf_name
    list_paf_path_against_ref = list(map(lambda name:
                            f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon/Round13.TGS_Gapcloser_Modified_V2_re.Unpl_Scf.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.Manual_Scaff.HiFi_trio.Hifiasm_Loc.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.ONT_trio.ONT_Amb.Manual_Scaff/Prealignment_for_Candidate_Filtering/Hap{hap}/Align/Align.REF_KOREFDraftASM_Hap{hap}_Split.QRY_{name}.asm5P.paf",
                            list_paf_name))

    list_paf_name_pairwise = list()
    list_paf_path_pairwise = list()

    path_paf_pairwise_format = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon/Round13.TGS_Gapcloser_Modified_V2_re.Unpl_Scf.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.Manual_Scaff.HiFi_trio.Hifiasm_Loc.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.ONT_trio.ONT_Amb.Manual_Scaff/Prealignment_for_Candidate_Filtering/Hap{hap}/Align/Align.REF.{name1}.QRY.{name2}.asm5P_triobin_{triobin}.paf"
    for ind1, name1 in enumerate(list_paf_name):
        for name2 in list(list_paf_name[ind1:]):
            n1, n2 = [name1, name2]
            path_paf = path_paf_pairwise_format.format(name1 = n1, name2 = n2, hap = hap, triobin = triobin)

            if not os.path.exists(path_paf):
                n1, n2 = [name2, name1]
                path_paf = path_paf_pairwise_format.format(name1 = n1, name2 = n2, hap = hap, triobin = triobin)
                assert os.path.exists(path_paf)

            list_paf_name_pairwise.append(f"{n1}::{n2}")
            list_paf_path_pairwise.append(path_paf)

    min_idy = 0.75
    min_align_contig_to_contig = 500
    min_align_contig_to_ref = 1000
    max_hang_contig_to_contig = 20000
    max_hang_contig_to_ref = 10000000
    min_align2hang_contig_to_contig = 0.75
    min_align2hang_contig_to_ref = 0.5
    n_alignedcontig_check = 0
    max_distance_ratio_from_last_contig = 0.5
    n_jobs = 40
# %%
