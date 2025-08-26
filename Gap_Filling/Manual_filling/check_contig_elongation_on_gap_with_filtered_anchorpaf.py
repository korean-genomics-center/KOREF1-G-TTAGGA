#%%
from pathlib import Path
from joblib import Parallel, delayed
from collections import Counter
import sys,os
from itertools import chain

import pandas as pd
import numpy as np

sys.path.append(str(Path(__file__).parents[0]))
from check_anchored_contig import read_all_paf_tables_with_filtering, get_candidate_alignment, filter_paf_of_end_aligned_query
from check_contig_overlap_elongation import read_all_pairwise_paf_tables_with_filtering, generate_contig_linkage_by_overlap, get_shortest_contig_paths_from_specific_contig_with_single_direction, get_statistics_of_elongated_path_on_paths_as_table, get_ranges_for_scaffold_elongation_from_path


def main(
    path_scaff_name,
    path_scaff_info,
    list_pafs_anchored,
    list_paf_name_pairwise,
    list_paf_path_pairwise,
    path_save_prefix,
    min_idy_contig_to_contig = 0.8,
    min_idy_contig_to_ref = 0.7,
    min_align_contig_to_contig = 500,
    min_align_contig_to_ref = 5000,
    max_hang_contig_to_contig = 10000,
    max_hang_contig_to_ref = 1000000,
    min_align2hang_contig_to_contig = 0.8,
    min_align2hang_contig_to_ref = 0.5,
    min_contig_length = 10000,
    n_alignedcontig_check = 50,
    max_distance_ratio_from_last_contig = 0.3,
    n_jobs = 20,
    max_elongation_contigs = 50
):
    os.makedirs(os.path.dirname(path_save_prefix), exist_ok=True)

    table_scaffname = pd.read_csv(path_scaff_name, sep = '\t', names = ["scaffname", "originalname"])
    dict_original2scaffname = dict(zip(table_scaffname["originalname"], table_scaffname["scaffname"].astype(str)))
    dict_scaff_info = read_scaff_info(path_scaff_info)
    list_contignames_between_gap = get_contig_pairs_with_gap(dict_scaff_info)

    table_all_refaligned_paf = pd.concat(list(map(lambda path: pd.read_csv(path, sep = '\t'), list_pafs_anchored)))
    table_all_refaligned_paf_filt = filter_paf_of_end_aligned_query(table_all_refaligned_paf, min_idy_contig_to_ref, min_align_contig_to_ref, max_hang_contig_to_ref, min_align2hang_contig_to_ref)
    table_all_refaligned_paf_filt["refname"] = table_all_refaligned_paf_filt["refname"].astype(str)
    table_all_refaligned_paf_filt = table_all_refaligned_paf_filt[table_all_refaligned_paf_filt["qrylen"] > min_contig_length]

    table_all_pairwise_paf_filt = read_all_pairwise_paf_tables_with_filtering(list_paf_name_pairwise, list_paf_path_pairwise, min_idy_contig_to_contig, min_align_contig_to_contig, max_hang_contig_to_contig, min_align2hang_contig_to_contig, n_jobs)
    table_all_pairwise_paf_filt = table_all_pairwise_paf_filt[np.logical_and(table_all_pairwise_paf_filt["reflen"] > min_contig_length, table_all_pairwise_paf_filt["qrylen"] > min_contig_length)]
    dict_contigs_redundant = get_similar_contig_sets_for_filter_out(table_all_pairwise_paf_filt)

    dict_contig_pairwise_alignment = generate_contig_linkage_by_overlap(table_all_pairwise_paf_filt)

    dict_chain_stats = find_elongation_path_from_aligned_contigs(list_contignames_between_gap, table_all_refaligned_paf_filt, dict_contig_pairwise_alignment, dict_contigs_redundant, n_alignedcontig_check, max_distance_ratio_from_last_contig, n_jobs, max_elongation_contigs)

    for contignames_key, table in dict_chain_stats.items():
        table.to_csv(
            f"{path_save_prefix}.contig_set_{contignames_key}.tsv",
            sep = '\t',
            index = False
        )

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

def get_similar_contig_sets_for_filter_out(table_all_pairwise_paf_filt):
    table_all_pairwise_paf_filt_highidy = table_all_pairwise_paf_filt[table_all_pairwise_paf_filt["identity"] > 0.95]
    table_all_pairwise_paf_filt_highidy = table_all_pairwise_paf_filt_highidy[table_all_pairwise_paf_filt_highidy["qryname"] != table_all_pairwise_paf_filt_highidy["refname"]]

    table_all_pairwise_paf_filt_highidy["ref_alignratio"] = table_all_pairwise_paf_filt_highidy.apply(lambda row: (row["ref_alignend"]-row["ref_alignstart"]) / row["reflen"], axis = 1)
    table_all_pairwise_paf_filt_highidy["qry_alignratio"] = table_all_pairwise_paf_filt_highidy.apply(lambda row: (row["qry_alignend"]-row["qry_alignstart"]) / row["qrylen"], axis = 1)

    table_all_pairwise_paf_filt_highidy_fulllength = table_all_pairwise_paf_filt_highidy[np.logical_and(table_all_pairwise_paf_filt_highidy["ref_alignratio"] > 0.95, table_all_pairwise_paf_filt_highidy["qry_alignratio"] > 0.95)]

    table_all_pairwise_paf_filt_highidy_fulllength["total_hang"] = table_all_pairwise_paf_filt_highidy_fulllength["qry_alignstart"] + table_all_pairwise_paf_filt_highidy_fulllength["ref_alignstart"] + table_all_pairwise_paf_filt_highidy_fulllength["qrylen"] - table_all_pairwise_paf_filt_highidy_fulllength["qry_alignend"] + table_all_pairwise_paf_filt_highidy_fulllength["reflen"] - table_all_pairwise_paf_filt_highidy_fulllength["ref_alignend"]
    table_all_pairwise_paf_filt_highidy_fulllength = table_all_pairwise_paf_filt_highidy_fulllength[table_all_pairwise_paf_filt_highidy_fulllength["total_hang"] < 1000]

    dict_contigs_redundant = dict()
    for _, row in table_all_pairwise_paf_filt_highidy_fulllength.iterrows():
        contig1 = row["qryname"]
        contig2 = row["refname"]

        if dict_contigs_redundant.get(contig1) == None: dict_contigs_redundant[contig1] = set()
        if dict_contigs_redundant.get(contig2) == None: dict_contigs_redundant[contig2] = set()

        dict_contigs_redundant[contig1].add(contig2)
        dict_contigs_redundant[contig2].add(contig1)
    return dict_contigs_redundant


def find_elongation_path_from_aligned_contigs(list_contignames_between_gap, table_all_refaligned_paf_filt, dict_contig_pairwise_alignment, dict_contigs_redundant, n_alignedcontig_check, max_distance_ratio_from_last_contig, n_jobs, max_elongation_contigs):
    dict_elongation_path = dict()

    for contignames in list_contignames_between_gap:
        contignames_key = "__to__".join(contignames)
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

        dict_contig_path_from_depart_to_dest = get_paths_from_depart_to_dest_contig_parallel(list_depart_contig_to_strand, dict_contig_dest_to_strand, dict_contig_pairwise_alignment, dict_contigs_redundant, n_jobs, max_elongation_contigs)
        # dict_contig_path_from_depart_to_dest = get_paths_from_depart_to_dest_contig_parallel(dict_contig_depart_to_strand, dict_contig_dest_to_strand, dict_contig_pairwise_alignment, n_jobs)
        table_path_stats = convert_depart_to_dest_paths_into_stats(dict_contig_path_from_depart_to_dest, dict_contig_pairwise_alignment, n_jobs)
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

def get_paths_from_depart_to_dest_contig_parallel(list_depart_contig_to_strand, dict_contig_dest_to_strand, dict_contig_alignment, dict_contigs_redundant, n_jobs, max_elongation_contigs):
    # list_contignames_depart = list(dict_contig_depart_to_strand.keys())
    list_contignames_dest = list(dict_contig_dest_to_strand.keys())

    with Parallel(n_jobs = n_jobs) as parallel:
        list_dict_paths_until_target = parallel(
            delayed(get_paths_until_another_aligned_contig_from_specific_contig)(
                contigname_depart__strand.split("__strand__")[0],
                list_contignames_dest,
                dict_contig_alignment,
                contigname_depart__strand.split("__strand__")[1],
                dict_contig_dest_to_strand,
                dict_contigs_redundant,
                max_elongation_contigs
            ) for contigname_depart__strand in list_depart_contig_to_strand
        )
    
    dict_path_between_depart_and_dest = dict()
    for contigname_depart__strand, dict_paths_until_dest in zip(list_depart_contig_to_strand, list_dict_paths_until_target):
        if len(dict_paths_until_dest.keys()) > 0:
            dict_path_between_depart_and_dest[contigname_depart__strand] = dict_paths_until_dest
    return dict_path_between_depart_and_dest

def get_paths_until_another_aligned_contig_from_specific_contig(contigname_start, list_contignames_end_target, dict_contig_alignment, contig_start_aligned_strand, dict_contig_dest_to_strand, dict_contigs_redundant, max_elongation_contigs):
    dict_path_from_start = get_shortest_contig_paths_from_specific_contig_with_single_direction(dict_contig_alignment, contigname_start, contig_start_aligned_strand == '-', max_dist = max_elongation_contigs, dict_contigs_redundant = dict_contigs_redundant)

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

            dict_paths_until_target[contig_end] = list_of_list_paths_with_correct_strand
    return dict_paths_until_target

def convert_depart_to_dest_paths_into_stats(dict_contig_path_from_depart_to_dest, dict_contig_pairwise_alignment, n_jobs):
    table_path_stats = pd.DataFrame()

    with Parallel(n_jobs = n_jobs) as parallel:
        list_of_list_tables = parallel(
            delayed(convert_specific_depart_contig_to_dest_path_into_stats)(
                contig_depart_strand.split("__strand__")[0],
                contig_depart_strand.split("__strand__")[1],
                dict_contig_path_from_depart_to_dest[contig_depart_strand],
                dict_contig_pairwise_alignment
            ) for contig_depart_strand in dict_contig_path_from_depart_to_dest.keys()
        )
    list_tables = list(chain(*list_of_list_tables))
    if len(list_tables) == 0:
        return pd.DataFrame()
    else:
        return pd.concat(list_tables)

def convert_specific_depart_contig_to_dest_path_into_stats(contig_depart, strand_depart, dict_path_to_dest, dict_contig_pairwise_alignment):
    strand_exchange = {
        '+' : '-',
        '-' : '+'
    }

    list_table_stats = list()
    for contig_dest in dict_path_to_dest.keys():
        list_paths_between_contig = dict_path_to_dest[contig_dest]
        if len(list_paths_between_contig) > 0:
            table_elongation_stat = get_statistics_of_elongated_path_on_paths_as_table(list_paths_between_contig, dict_contig_pairwise_alignment, "end" if strand_depart == '+' else "front")
            table_elongation_stat["Depart"] = contig_depart
            table_elongation_stat["Destination"] = contig_dest

            list_dict_contig_path_trimming = table_elongation_stat[["Path", "Direction"]].apply(
                lambda row: get_ranges_for_scaffold_elongation_from_path(row["Path"].split("->"), dict_contig_pairwise_alignment, row["Direction"]),
                axis = 1
            )

            table_elongation_stat["Length"] = list(map(lambda dict_contig_trimming: get_total_path_length(dict_contig_trimming), list_dict_contig_path_trimming))
            table_elongation_stat["Depart_strand"] = strand_depart

            if strand_depart == "+":
                table_elongation_stat["Destination_strand"] = list(map(lambda dict_contig_trimming: dict_contig_trimming[contig_dest]["strand"], list_dict_contig_path_trimming))
            else:
                table_elongation_stat["Destination_strand"] = list(map(lambda dict_contig_trimming: strand_exchange[dict_contig_trimming[contig_dest]["strand"]], list_dict_contig_path_trimming))

            list_table_stats.append(table_elongation_stat)
    return list_table_stats



def get_total_path_length(dict_ranges_for_trimming):
    list_all_ranges = list(map(lambda trimming_config: trimming_config["range"], dict_ranges_for_trimming.values()))
    list_range_length = list(map(lambda ranges: ranges[1]-ranges[0]+1, list_all_ranges))

    return sum(list_range_length)


#%%
if __name__ == "__main__":
    list_paf_name = [
        "Hifiasm_Binned_HiFi",
        "Hifiasm_Binned_ALT_HiFi",
        "Hifiasm_Binned_HiFi_with0",
        "Hifiasm_Binned_ALT_HiFi_with0",
        "Flye_Binned_HiFi_with0",
        "Flye_Binned_HiFi",
        "Unplaced_Scaffold"
    ]

    hap = '1'
    triobin = 'p'
    
    path_scaff_name = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon_wo_RawONT/Round8.Manual_Scaff/Hap{hap}/Prealignment/tgsseqsplit/CuratedDraftASM.KOREF_S1.Trio_Specific.hap{hap}.Round8.GapFilling.Manual_Scaff.LocalASM.Unplaced_Scf.Hifiasm_Hifi_Bin.Hifiasm_Hifi_Bin0.Flye_Hifi_Bin.Flye_Hifi_Bin0.Manual_Scaff.name_map"
    path_scaff_info = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon_wo_RawONT/Round8.Manual_Scaff/Hap{hap}/Prealignment/tgsseqsplit/CuratedDraftASM.KOREF_S1.Trio_Specific.hap{hap}.Round8.GapFilling.Manual_Scaff.LocalASM.Unplaced_Scf.Hifiasm_Hifi_Bin.Hifiasm_Hifi_Bin0.Flye_Hifi_Bin.Flye_Hifi_Bin0.Manual_Scaff.orignial_scaff_infos"

    list_pafs_anchored = [
        "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon_wo_RawONT/Round8.Manual_Scaff/Hap1/Prealignment/Anchored_Alignments/Anchored_Alignments.REF_KOREFDraftASM_Hap1_Split.chrname_chr22.contig_index_1.End.tsv",
        "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon_wo_RawONT/Round8.Manual_Scaff/Hap1/Prealignment/Anchored_Alignments/Anchored_Alignments.REF_KOREFDraftASM_Hap1_Split.chrname_chr22.contig_index_2.Front.tsv",
        "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon_wo_RawONT/Round8.Manual_Scaff/Hap1/Prealignment/Anchored_Alignments/Anchored_Alignments.REF_KOREFDraftASM_Hap1_Split.chrname_chr2.contig_index_1.End.tsv",
        "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon_wo_RawONT/Round8.Manual_Scaff/Hap1/Prealignment/Anchored_Alignments/Anchored_Alignments.REF_KOREFDraftASM_Hap1_Split.chrname_chr2.contig_index_2.Front.tsv",
        "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon_wo_RawONT/Round8.Manual_Scaff/Hap1/Prealignment/Anchored_Alignments/Anchored_Alignments.REF_KOREFDraftASM_Hap1_Split.chrname_chr7.contig_index_1.End.tsv",
        "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon_wo_RawONT/Round8.Manual_Scaff/Hap1/Prealignment/Anchored_Alignments/Anchored_Alignments.REF_KOREFDraftASM_Hap1_Split.chrname_chr7.contig_index_2.Front.tsv",
        "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon_wo_RawONT/Round8.Manual_Scaff/Hap1/Prealignment/Anchored_Alignments/Anchored_Alignments.REF_KOREFDraftASM_Hap1_Split.chrname_chr8.contig_index_1.End.tsv",
        "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon_wo_RawONT/Round8.Manual_Scaff/Hap1/Prealignment/Anchored_Alignments/Anchored_Alignments.REF_KOREFDraftASM_Hap1_Split.chrname_chr8.contig_index_2.Front.tsv"
    ]

    list_paf_name_pairwise = list()
    list_paf_path_pairwise = list()

    path_paf_pairwise_format = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon_wo_RawONT/Contig_Alignment/Hap{hap}/Align.Pairwise/Align.Hap{hap}.REF.{name1}.QRY.{name2}.asm5Pc_triobin_{triobin}.filt_match100.paf"
    for name1 in list_paf_name:
        for name2 in list_paf_name:
            n1, n2 = [name1, name2]
            path_paf = path_paf_pairwise_format.format(name1 = n1, name2 = n2, hap = hap, triobin = triobin)

            if not os.path.exists(path_paf):
                n1, n2 = [name2, name1]
                path_paf = path_paf_pairwise_format.format(name1 = n1, name2 = n2, hap = hap, triobin = triobin)
                assert os.path.exists(path_paf)

            list_paf_name_pairwise.append(f"{n1}::{n2}")
            list_paf_path_pairwise.append(path_paf)

    path_save_prefix = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon_wo_RawONT/Round8.Manual_Scaff/Hap{hap}/Prealignment/Align_Graph/Alignment_Graph.stats_of_contigs_joining_gap"

    min_idy_contig_to_contig = 0.9
    min_idy_contig_to_ref = 0.9
    min_align_contig_to_contig = 500
    min_align_contig_to_ref = 1000
    max_hang_contig_to_contig = 10000
    max_hang_contig_to_ref = 10000
    min_align2hang_contig_to_contig = 0.9
    min_align2hang_contig_to_ref = 0.9
    min_contig_length = 0
    n_alignedcontig_check = 0
    max_distance_ratio_from_last_contig = 0.3
    n_jobs = 40
    max_elongation_contigs = 1000

    # main(
    #     path_scaff_name,
    #     path_scaff_info,
    #     list_pafs_anchored,
    #     list_paf_name_pairwise,
    #     list_paf_path_pairwise,
    #     path_save_prefix,
    #     min_idy_contig_to_contig = min_idy_contig_to_contig,
    #     min_idy_contig_to_ref = min_idy_contig_to_ref,
    #     min_align_contig_to_contig = min_align_contig_to_contig,
    #     min_align_contig_to_ref = min_align_contig_to_ref,
    #     max_hang_contig_to_contig = max_hang_contig_to_contig,
    #     max_hang_contig_to_ref = max_hang_contig_to_ref,
    #     min_align2hang_contig_to_contig = min_align2hang_contig_to_contig,
    #     min_align2hang_contig_to_ref = min_align2hang_contig_to_ref,
    #     min_contig_length = min_contig_length,
    #     n_alignedcontig_check = n_alignedcontig_check,
    #     max_distance_ratio_from_last_contig = max_distance_ratio_from_last_contig,
    #     n_jobs = n_jobs,
    #     max_elongation_contigs = max_elongation_contigs
    # )
# %%
