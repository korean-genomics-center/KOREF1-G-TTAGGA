#%%
import os,sys
from copy import deepcopy
import heapq

import pandas as pd
import numpy as np


class ContigNode(object):
    def __init__(self, contigname, contiglength):
        self.contigname = contigname
        self.contiglength = contiglength
        self.front_overlap = dict()
        self.end_overlap = dict()
    
    def get_next_node(self, get_front = False):
        if get_front:
            return self.front_overlap
        else:
            return self.end_overlap

    def add_front_alignment(self, other_contig_name, align_info, this_hang, other_hang):
        if self.front_overlap.get(other_contig_name) != None and self.front_overlap[other_contig_name]["info"].is_this_better(align_info):
            return
        self.front_overlap[other_contig_name] = {
            "info":align_info,
            "hang":{
                "this":this_hang,
                "other":other_hang
            }
        }

    def add_end_alignment(self, other_contig_name, align_info, this_hang, other_hang):
        if self.end_overlap.get(other_contig_name) != None and self.end_overlap[other_contig_name]["info"].is_this_better(align_info):
            return
        self.end_overlap[other_contig_name] = {
            "info":align_info,
            "hang":{
                "this":this_hang,
                "other":other_hang
            }
        }
    
    def __repr__(self):
        text = f"Contig: {self.contigname}\n"
        text += f"Length: {self.contiglength}\n"
        text += "Front: \n"
        for key, val in self.front_overlap.items():
            text += f"\tContig: {key}\n"
            text += f"\t\t{val['info']}\n"
            text += f"\t\tHang: {val['hang']}\n"
        text += "\n"
        text += "End: \n"
        for key, val in self.end_overlap.items():
            text += f"\tContig: {key}\n"
            text += f"\t\t{val['info']}\n"
            text += f"\t\tHang: {val['hang']}\n"
        return text
    
    def get_dictionarized_object(self):
        dict_this = {
            "Front":dict(),
            "End":dict()
        }
        
        for key, val in self.front_overlap.items():
            dict_this["Front"][key] = {
                "Alignment" : val["info"].__dict__,
                "Hang" : val["hang"]
            }
        for key, val in self.end_overlap.items():
            dict_this["End"][key] = {
                "Alignment" : val["info"].__dict__,
                "Hang" : val["hang"]
            }
        return dict_this


class AlignmentINFO(object):
    def __init__(self, lenalign, idy, alignratio, alignstrand, refname, qryname, reflen, qrylen, refalign, qryalign):
        self.lenalign = lenalign
        self.identity = idy
        self.alignratio = alignratio
        self.alignstrand = alignstrand
        self.refname = refname
        self.qryname = qryname
        self.reflen = reflen
        self.qrylen = qrylen
        self.refalign = refalign
        self.qryalign = qryalign
    
    def is_this_better(self, other):
        if self.identity > 0.95 and other.identity > 0.95:
            if self.lenalign > other.lenalign:
                return True
        elif self.identity > other.identity:
                return True
        return False

    def __repr__(self):
        return f"Strand: {self.alignstrand}, Alignment Length: {self.lenalign}, Identity: {self.identity}, Alignment Ratio: {self.alignratio}"


def run(
    list_paf_namepair,
    list_paf_path,
    min_idy,
    min_align,
    max_hang,
    min_alignratio
):
    table_all_paf_filt = read_all_pairwise_paf_tables_with_filtering(list_paf_namepair, list_paf_path, min_idy, min_align, max_hang, min_alignratio)
    dict_contig_alignment = generate_contig_linkage_by_overlap(table_all_paf_filt)
    dict_contig_to_elongated_paths = get_all_elongation_paths_from_all_contig(dict_contig_alignment)
    dict_contig_alignment_for_json = dict()
    for contigname, contigobj in dict_contig_alignment.items():
        dict_contig_alignment_for_json[contigname] = contigobj.get_dictionarized_object()

def read_all_pairwise_paf_tables_with_filtering(list_paf_namepair, list_paf_path, min_idy, min_align, max_hang, min_align2hang):
    list_table_paf = list()
    for namepair, path_paf in zip(list_paf_namepair, list_paf_path):
        assert os.path.exists(path_paf), path_paf
        table_paf_name = read_paf(path_paf)
        table_paf_name_filt_end_aligned = filter_paf_of_end_aligned_query(table_paf_name, min_idy, min_align, max_hang, min_align2hang)
        name_ref, name_qry = namepair.split('::')
        table_paf_name_filt_end_aligned.loc[:, "qryname"] = table_paf_name_filt_end_aligned["qryname"].apply(lambda val: f"{name_qry}::{val}")
        table_paf_name_filt_end_aligned.loc[:, "refname"] = table_paf_name_filt_end_aligned["refname"].apply(lambda val: f"{name_ref}::{val}")
        list_table_paf.append(table_paf_name_filt_end_aligned)
    return pd.concat(list_table_paf)

def read_paf(path_paf):
    table_paf = pd.read_csv(path_paf, sep = '\t', names = list(range(1, 13)), usecols = list(range(1, 13)))
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

def filter_paf_of_end_aligned_query(table_paf, min_idy, min_align, max_hang, min_align2hang):
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

def generate_contig_linkage_by_overlap(table_paf):
    dict_contigs = dict()
    for _, row in table_paf.iterrows():
        dict_alignment_info = row.to_dict()

        align_strand = dict_alignment_info["strand"]
        contig1_name = dict_alignment_info["refname"]
        contig2_name = dict_alignment_info["qryname"]

        if dict_contigs.get(contig1_name) == None:
            dict_contigs[contig1_name] = ContigNode(contig1_name, dict_alignment_info["reflen"])
        if dict_contigs.get(contig2_name) == None:
            dict_contigs[contig2_name] = ContigNode(contig2_name, dict_alignment_info["qrylen"])
        if contig1_name == contig2_name: continue

        contig1_hang, contig2_hang, is_contig1_left = check_align_hang(dict_alignment_info, align_strand == '-')
        if contig1_hang == None:
            continue
        
        align_info = AlignmentINFO(dict_alignment_info["align"], dict_alignment_info["identity"], dict_alignment_info["aligned2hang"], align_strand, contig1_name, contig2_name, dict_alignment_info["reflen"], dict_alignment_info["qrylen"], (dict_alignment_info["ref_alignstart"], dict_alignment_info["ref_alignend"]), (dict_alignment_info["qry_alignstart"], dict_alignment_info["qry_alignend"]))
        if is_contig1_left:
            dict_contigs[contig1_name].add_end_alignment(contig2_name, align_info, contig1_hang, contig2_hang)
            if align_strand == '+':
                dict_contigs[contig2_name].add_front_alignment(contig1_name, align_info, contig2_hang, contig1_hang)
            else:
                dict_contigs[contig2_name].add_end_alignment(contig1_name, align_info, contig2_hang, contig1_hang)
        else:
            dict_contigs[contig1_name].add_front_alignment(contig2_name, align_info, contig1_hang, contig2_hang)
            if align_strand == '+':
                dict_contigs[contig2_name].add_end_alignment(contig1_name, align_info, contig2_hang, contig1_hang)
            else:
                dict_contigs[contig2_name].add_front_alignment(contig1_name, align_info, contig2_hang, contig1_hang)
    return dict_contigs

def check_align_hang(dict_alignment_info, is_rev):
    contig1_lefthang = dict_alignment_info["ref_alignstart"]
    contig1_righthang = dict_alignment_info["reflen"] - dict_alignment_info["ref_alignend"]

    if is_rev:
        contig2_lefthang = dict_alignment_info["qrylen"] - dict_alignment_info["qry_alignend"]
        contig2_righthang = dict_alignment_info["qry_alignstart"]
    else:
        contig2_lefthang = dict_alignment_info["qry_alignstart"]
        contig2_righthang = dict_alignment_info["qrylen"] - dict_alignment_info["qry_alignend"]
    
    contig1_hang = None
    contig2_hang = None
    is_contig1_left = True
    if contig1_lefthang > contig2_lefthang and contig1_righthang < contig2_righthang:
        contig1_hang = contig1_lefthang
        contig2_hang = contig2_righthang
    elif contig1_lefthang < contig2_lefthang and contig1_righthang > contig2_righthang:
        contig1_hang = contig1_righthang
        contig2_hang = contig2_lefthang
        is_contig1_left = False
    
    return (contig1_hang, contig2_hang, is_contig1_left)

def get_all_elongation_paths_from_all_contig(dict_contig_align):
    dict_contig_to_elongation_paths = dict()
    for contig_start in dict_contig_align.keys():
        front_paths_to_other_contigs = get_shortest_contig_paths_from_specific_contig_with_single_direction(dict_contig_align, contig_start, True)
        end_paths_to_other_contigs = get_shortest_contig_paths_from_specific_contig_with_single_direction(dict_contig_align, contig_start, False)

        front_uniquely_elongated_paths = remain_maximum_unique_elongated_path_only(front_paths_to_other_contigs)
        end_uniquely_elongated_paths = remain_maximum_unique_elongated_path_only(end_paths_to_other_contigs)

        dict_contig_to_elongation_paths[contig_start] = {
            "Elongate_Front" : front_uniquely_elongated_paths,
            "Elongate_End" : end_uniquely_elongated_paths
        }
    return dict_contig_to_elongation_paths

def get_shortest_contig_paths_from_specific_contig_with_single_direction(dict_contig_align, start_contig, go_front = True, max_dist = 20, dict_contigs_redundant = None):
    graph_queue = list()
    heapq.heappush(graph_queue, (0, start_contig, [start_contig], go_front))
    shortest_paths = dict()
    shortest_distance = dict()

    shortest_distance[start_contig] = 0

    while graph_queue:
        dist, curr_contig, path, check_front = heapq.heappop(graph_queue)

        if curr_contig in shortest_distance and dist > shortest_distance[curr_contig]:
            continue

        if curr_contig not in shortest_paths:
            shortest_paths[curr_contig] = list()
        shortest_paths[curr_contig].append(path)

        if dict_contig_align.get(curr_contig) == None:
            print(f"No alignment in {curr_contig}", flush = True)
            continue
        
        dict_next_contigs = dict_contig_align[curr_contig].get_next_node(get_front = check_front)

        list_next_contigs = list(dict_next_contigs.keys())
        if dict_contigs_redundant != None:
            dict_redundancy_of_next_contigs = get_contigs_grouped_by_redundancy(list_next_contigs, dict_contigs_redundant)
            list_next_contigs = get_best_identity_from_redundancy_group(dict_redundancy_of_next_contigs, dict_next_contigs)

        for neighbor_contig in list_next_contigs:
            dist_new = dist + 1

            if dist_new > max_dist:
                continue

            check_front_of_neighbor = ((int(check_front)-0.5) * (1 if dict_next_contigs[neighbor_contig]["info"].alignstrand == '+' else -1)) > 0
            if neighbor_contig not in shortest_distance or dist_new < shortest_distance[neighbor_contig]:
                shortest_distance[neighbor_contig] = dist_new
                heapq.heappush(graph_queue, (dist_new, neighbor_contig, path+[neighbor_contig], check_front_of_neighbor))
            elif dist_new == shortest_distance[neighbor_contig]:
                heapq.heappush(graph_queue, (dist_new, neighbor_contig, path+[neighbor_contig], check_front_of_neighbor))
    
    return shortest_paths

def get_contigs_grouped_by_redundancy(list_contigs, dict_contigs_redundant):
    dict_redundancy_group = dict()

    curr_group = 0
    set_contigs_left = set(list_contigs)
    while len(set_contigs_left) > 0:
        contigname = set_contigs_left.pop()
        if dict_contigs_redundant.get(contigname) != None:
            other_contigs_in_redundancy = dict_contigs_redundant[contigname] & set_contigs_left
            dict_redundancy_group[curr_group] = [contigname] + list(other_contigs_in_redundancy)

            set_contigs_left -= other_contigs_in_redundancy
        else:
            dict_redundancy_group[curr_group] = [contigname]
        curr_group += 1
    return dict_redundancy_group

def get_best_identity_from_redundancy_group(dict_redundancy_of_contigs, dict_contig_align):
    list_contigs_selected = list()

    for _, list_contigs_in_redundancy in dict_redundancy_of_contigs.items():
        list_identity_of_contigs = list(map(lambda contigname: dict_contig_align[contigname]["info"].identity, list_contigs_in_redundancy))
        dict_contig_to_idy = dict(zip(list_contigs_in_redundancy, list_identity_of_contigs))
        sorted_contignames_in_redundancy_by_idy = sorted(list_contigs_in_redundancy, key = dict_contig_to_idy.__getitem__, reverse = True)

        list_contigs_selected.append(sorted_contignames_in_redundancy_by_idy[0])
    return list_contigs_selected

def remain_maximum_unique_elongated_path_only(dict_shortest_paths):
    list_all_paths = list()
    for _, list_paths in dict_shortest_paths.items():
        for paths in list_paths:
            list_all_paths.append(paths)

    list_maximum_unique_paths = list()
    for path in list_all_paths:
        is_path_uniquely_elongated = len(path) > 1
        for path_check in list_all_paths:
            if len(path) == len(path_check):
                continue
            if check_list_startswith_this_list(path_check, path):
                is_path_uniquely_elongated = False
        if is_path_uniquely_elongated:
            list_maximum_unique_paths.append(path)
    return list_maximum_unique_paths
                
def check_list_startswith_this_list(list_check, list_qry):
    if len(list_check) < len(list_qry):
        return False
    for ind, val in enumerate(list_qry):
        if list_check[ind] != val:
            return False
    return True

def get_statistics_of_all_elongation_path_as_table(contig_start, dict_contig_to_elongated_paths, dict_contig_alignment):
    table_elongation_stat = pd.DataFrame(
        columns = ["Path", "Direction", "Length", "Mean_Identity", "Min_Identity", "Max_Identity", "Mean_AlignmentRatio", "Min_AlignmentRatio", "Max_AlignmentRatio"]
    )
    
    for direction in ["Front", "End"]:
        table_elongation_stat_of_direction = get_statistics_of_elongated_path_on_paths_as_table(
            dict_contig_to_elongated_paths[contig_start][f"Elongate_{direction}"],
            dict_contig_alignment,
            direction
        )
        table_elongation_stat = pd.concat([table_elongation_stat, table_elongation_stat_of_direction])

    return table_elongation_stat

def get_statistics_of_elongated_path_on_paths_as_table(list_paths, dict_contig_alignment, direction):
    table_elongation_stat = pd.DataFrame(
        columns = ["Path", "Direction", "Length", "Mean_Identity", "Min_Identity", "Max_Identity", "Mean_AlignmentRatio", "Min_AlignmentRatio", "Max_AlignmentRatio"]
    )

    for path in list_paths:
        dict_stats_of_path = get_statistics_of_elongated_path(path, dict_contig_alignment, direction.lower())

        series_stats = pd.Series({
            "Path":'->'.join(path), 
            "Direction":direction, 
            "Length":sum(dict_stats_of_path["length"]), 
            "Mean_Identity":apply_function_with_empty_list_error_case(dict_stats_of_path["identity"], np.mean), 
            "Min_Identity":apply_function_with_empty_list_error_case(dict_stats_of_path["identity"], min), 
            "Max_Identity":apply_function_with_empty_list_error_case(dict_stats_of_path["identity"], max), 
            "Mean_AlignmentRatio":apply_function_with_empty_list_error_case(dict_stats_of_path["alignratio"], np.mean), 
            "Min_AlignmentRatio":apply_function_with_empty_list_error_case(dict_stats_of_path["alignratio"], min), 
            "Max_AlignmentRatio":apply_function_with_empty_list_error_case(dict_stats_of_path["alignratio"], max)
        })
        table_elongation_stat = table_elongation_stat.append(series_stats, ignore_index = True)
    return table_elongation_stat

def apply_function_with_empty_list_error_case(values, func):
    if len(values) == 0:
        return np.nan
    else:
        return func(values)

def get_statistics_of_elongated_path(list_elongated_contignames, dict_contig_alignment, direction):

    list_elongated_length = list()
    list_identity = list()
    list_alignment_ratio = list()

    direction_exchange = {
        "front":"end",
        "end":"front"
    }

    current_direction = direction
    current_contig = list_elongated_contignames[0]
    for contig_elong in list_elongated_contignames[1:]:
        elong_alignment = dict_contig_alignment[current_contig].__dict__[f"{current_direction}_overlap"][contig_elong]

        list_elongated_length.append(elong_alignment["hang"]["other"])
        list_identity.append(elong_alignment["info"].identity)
        list_alignment_ratio.append(elong_alignment["info"].alignratio)

        if elong_alignment["info"].alignstrand == '-':
            current_direction = direction_exchange[current_direction]
        current_contig = contig_elong

    dict_stats = {
        "length":list_elongated_length,
        "identity":list_identity,
        "alignratio":list_alignment_ratio
    }

    return dict_stats

def get_ranges_for_scaffold_elongation_from_path(list_elongated_contignames, dict_contig_alignment, direction):
    dict_contig_to_ranges = dict()

    strand_exchange = {
        '+' : '-',
        '-' : '+'
    }
    direction_exchange = {
        "front":"end",
        "end":"front"
    }
    
    current_direction = direction
    current_strand = '+'
    current_contig = list_elongated_contignames[0]

    current_contig_range = [1, dict_contig_alignment[current_contig].contiglength]

    for elong_contig in list_elongated_contignames[1:]:
        elong_alignment = dict_contig_alignment[current_contig].__dict__[f"{current_direction}_overlap"][elong_contig]

        assert elong_alignment["info"].refname == current_contig or elong_alignment["info"].qryname == current_contig
        assert elong_alignment["info"].refname == elong_contig or elong_alignment["info"].qryname == elong_contig
        aligninfo_currcontig_dataset = "ref" if elong_alignment["info"].refname == current_contig else "qry"
        aligninfo_elongcontig_dataset = "ref" if elong_alignment["info"].refname == elong_contig else "qry"

        current_contig_align = elong_alignment["info"].__dict__[f"{aligninfo_currcontig_dataset}align"]
        elong_contig_align = elong_alignment["info"].__dict__[f"{aligninfo_elongcontig_dataset}align"]

        if current_strand == '+':
            if direction == "front":
                current_contig_range[0] = current_contig_align[0]
            else:
                current_contig_range[1] = current_contig_align[1]
        else:
            if direction == "front":
                current_contig_range[1] = current_contig_align[1]
            else:
                current_contig_range[0] = current_contig_align[0]

        dict_contig_to_ranges[current_contig] = {
            "strand":current_strand,
            "range":current_contig_range
        }

        if elong_alignment["info"].alignstrand == '-':
            current_direction = direction_exchange[current_direction]
            current_strand = strand_exchange[current_strand]
            
        current_contig_range = [1, dict_contig_alignment[elong_contig].contiglength]
        if current_strand == '+':
            if direction == "front":
                current_contig_range[1] = elong_contig_align[0]
            else:
                current_contig_range[0] = elong_contig_align[1]
        else:
            if direction == "front":
                current_contig_range[0] = elong_contig_align[1]
            else:
                current_contig_range[1] = elong_contig_align[0]
        current_contig = elong_contig
    
    dict_contig_to_ranges[current_contig] = {
        "strand":current_strand,
        "range" : current_contig_range
    }
    return dict_contig_to_ranges


# %%
if __name__ == "__main__":
    hap = 2
    triobin = 'm'

    list_names_contig = [
        "Hifiasm_Binned_pri",
        "Flye_HiFi_Binned",
        "Flye_ONT_Binned",
        "Unplaced_Scaffold"
    ]

    list_paf_namepair = list()
    list_paf_path = list()

    for ind1, name1 in enumerate(list_names_contig):
        for name2 in list(list_names_contig[ind1:]):
            path_paf = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon/Round13.TGS_Gapcloser_Modified_V2_re.Unpl_Scf.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.Manual_Scaff.HiFi_trio.Hifiasm_Loc.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.ONT_trio.ONT_Amb.Manual_Scaff/Prealignment_for_Candidate_Filtering/Hap{hap}/Align/Align.REF.{name1}.QRY.{name2}.asm5P_triobin_{triobin}.paf"

            list_paf_namepair.append(f"{name1}::{name2}")
            list_paf_path.append(path_paf)

    min_idy = 0.8
    min_align = 500
    max_hang = 10000
    min_alignratio = 0.8

    # run(
    #     list_paf_namepair,
    #     list_paf_path,
    #     min_idy,
    #     min_align,
    #     max_hang,
    #     min_alignratio
    # )
# %%
