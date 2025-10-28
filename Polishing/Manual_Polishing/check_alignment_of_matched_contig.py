#%%
import os,sys,re,math
from copy import deepcopy
from joblib import Parallel, delayed
from collections import Counter
from itertools import chain

import pandas as pd
import numpy as np


def run(
    path_paf,
    path_lowqv_region_info,
    path_output,
    anchor_min_idy,
    anchor_min_align,
    anchor_max_wronghang,
    anchor_min_align2hang,
    n_jobs
):
    dict_lowqv_region = read_lowqv_anchor_region(path_lowqv_region_info)
    dict_qry_to_anchor_aligned = read_paf_considering_anchor(path_paf, dict_lowqv_region)
    dict_qry_to_anchor_align_table = collect_alignments_on_anchor_with_filtering_parallel(dict_qry_to_anchor_aligned, dict_lowqv_region, anchor_min_idy, anchor_min_align, anchor_max_wronghang, anchor_min_align2hang, n_jobs)

    table_all_anchor_align = pd.concat(dict_qry_to_anchor_align_table.values())
    table_all_anchor_align.to_csv(path_output, sep = '\t', index = False)

def read_lowqv_anchor_region(path_lowqv_region_info):
    table_lowqv_region = pd.read_csv(path_lowqv_region_info, sep = '\t')
    table_lowqv_region["name"] = table_lowqv_region.apply(lambda row: f"::{row['chr']}:{row['range_start']}-{row['range_end']}", axis = 1)

    assert len(table_lowqv_region["name"].unique()) == table_lowqv_region.shape[0]

    table_lowqv_region["front_anchor_start"] = 0
    table_lowqv_region["front_anchor_end"] = table_lowqv_region["lowqv_start"] - table_lowqv_region["range_start"] - 1
    table_lowqv_region["tail_anchor_start"] = table_lowqv_region["lowqv_end"] - table_lowqv_region["range_start"]
    table_lowqv_region["tail_anchor_end"] = table_lowqv_region["range_end"] - table_lowqv_region["range_start"] - 1

    dict_name_to_lowqv_info = dict()
    for _, row in table_lowqv_region.iterrows():
        dict_name_to_lowqv_info[row["name"]] = {
            "anchor_front": {
                "start" : row["front_anchor_start"],
                "end" : row["front_anchor_end"]
            },
            "anchor_tail" : {
                "start" : row["tail_anchor_start"],
                "end" : row["tail_anchor_end"]
            },
            "lowqv_region" : {
                "start" : row["front_anchor_end"] + 1,
                "end" : row["tail_anchor_start"] - 1
            }
        }
    
    return dict_name_to_lowqv_info

def read_paf_considering_anchor(path_paf, dict_lowqv_region):
    dict_qry_to_anchor_aligned = dict()

    with open(path_paf, 'r') as fr:
        for line in fr:
            list_values = line.strip('\n').split('\t')

            dict_required_values = dict(zip(
                ["qryname", "qrylen", "qry_alignstart", "qry_alignend", "strand", "refname", "reflen", "ref_alignstart", "ref_alignend", "matchlen", "alignlen"],
                list_values[:11]
            ))

            qry_alignstart = int(dict_required_values["qry_alignstart"])
            qry_alignend = int(dict_required_values["qry_alignend"])-1

            is_front_anchor_included = check_is_align_include_anchor(qry_alignstart, qry_alignend, dict_lowqv_region[dict_required_values["qryname"]]["anchor_front"]["start"], dict_lowqv_region[dict_required_values["qryname"]]["anchor_front"]["end"])
            is_tail_anchor_included = check_is_align_include_anchor(qry_alignstart, qry_alignend, dict_lowqv_region[dict_required_values["qryname"]]["anchor_tail"]["start"], dict_lowqv_region[dict_required_values["qryname"]]["anchor_tail"]["end"])

            if not (is_front_anchor_included or is_tail_anchor_included):
                continue

            if dict_qry_to_anchor_aligned.get(dict_required_values["qryname"]) == None:
                dict_qry_to_anchor_aligned[dict_required_values["qryname"]] = dict()

            if dict_qry_to_anchor_aligned[dict_required_values["qryname"]].get(dict_required_values["refname"]) == None:
                dict_qry_to_anchor_aligned[dict_required_values["qryname"]][dict_required_values["refname"]] = list()
            
            dict_anchor_including = {
                "front" : is_front_anchor_included,
                "tail" : is_tail_anchor_included
            }

            val_cigar = list(filter(lambda val: val.startswith("cg:Z:"), list_values))
            assert len(val_cigar) == 1
            val_cigar = val_cigar[0].replace("cg:Z:", '')

            dict_qry_to_anchor_aligned[dict_required_values["qryname"]][dict_required_values["refname"]].append(
                {
                    "anchor_including" : dict_anchor_including,
                    "alignment" : dict_required_values,
                    "cigar" : val_cigar
                }
            )
    return dict_qry_to_anchor_aligned

def check_is_align_include_anchor(alignstart, alignend, anchor_start, anchor_end):
    if alignstart >= anchor_start:
        if alignstart < anchor_end:
            return True
        else:
            return False
    elif alignend > anchor_start:
        return True
    else:
        return False
    
def collect_alignments_on_anchor_with_filtering_parallel(dict_qry_to_anchor_aligned, dict_lowqv_region, anchor_min_idy, anchor_min_align, anchor_max_wronghang, anchor_min_align2hang, n_jobs):
    list_qrys = list(dict_qry_to_anchor_aligned.keys())
    
    # if n_jobs < 2:
    #     list_tables_by_qry = []
    #     for qry in list_qrys:
    #         result = get_alignment_of_single_interest_anchor(
    #             dict_qry_to_anchor_aligned[qry],
    #             dict_lowqv_region[qry],
    #             anchor_min_align, anchor_min_align2hang, anchor_max_wronghang, anchor_min_idy
    #         )
    #         list_tables_by_qry.append(result)
    # else:
    #     with Parallel(n_jobs = n_jobs) as parallel:
    #         list_tables_by_qry = parallel(
    #             delayed(get_alignment_of_single_interest_anchor)(
    #                 dict_qry_to_anchor_aligned[qry],
    #                 dict_lowqv_region[qry],
    #                 anchor_min_align, anchor_min_align2hang, anchor_max_wronghang, anchor_min_idy
    #             ) for qry in list_qrys
    #         )
    
    list_tables_by_qry = []
    for qry in list_qrys:
        result = get_alignment_of_single_interest_anchor(
            dict_qry_to_anchor_aligned[qry],
            dict_lowqv_region[qry],
            anchor_min_align, anchor_min_align2hang, anchor_max_wronghang, anchor_min_idy, n_jobs
        )
        list_tables_by_qry.append(result)
    
    dict_qry_to_anchor_align_table = dict(zip(list_qrys, list_tables_by_qry))
    return dict_qry_to_anchor_align_table

def get_alignment_of_single_interest_anchor(dict_anchor_align, dict_region, anchor_min_align, anchor_min_align2hang, anchor_max_wronghang, anchor_min_idy, n_jobs):
    list_all_alignments = list(chain(*dict_anchor_align.values()))
    
    list_splitted_all_alignments = split_list_into_n_list_values(list_all_alignments, min(n_jobs, len(list_all_alignments)))
    
    with Parallel(n_jobs = n_jobs) as parallel:
        list_filtered_anchor_paf = parallel(delayed(make_table_and_filter_of_list_alignments_on_anchor)(
            list_alignments_batch,
            dict_region, anchor_min_align, anchor_min_align2hang, anchor_max_wronghang, anchor_min_idy
        )for list_alignments_batch in list_splitted_all_alignments)
    table_anchor_paf_filtered = pd.concat(list_filtered_anchor_paf)
    
    # table_anchor_paf = make_table_of_anchor_including_alignments(dict_anchor_align)
    # table_anchor_paf_filtered = filter_table_with_anchor_alignment_position(table_anchor_paf, dict_region, anchor_min_align, anchor_min_align2hang, anchor_max_wronghang, anchor_min_idy)
    return table_anchor_paf_filtered

def split_list_into_n_list_values(list_values, n_split):
    n_values = len(list_values)
    n_values_per_split = n_values // n_split
    leftover_values = n_values % n_split
    
    list_split_values = list()
    ind_start = 0
    for i in range(n_split):
        n_values_thisind = n_values_per_split + int(i < leftover_values)
        list_split_values.append(list_values[ind_start:ind_start+n_values_thisind])
        
        ind_start += n_values_thisind
    return list_split_values

def make_table_and_filter_of_list_alignments_on_anchor(list_dict_alignments, dict_region, anchor_min_align, anchor_min_align2hang, anchor_max_wronghang, anchor_min_idy):
    table_anchor_paf = make_table_of_anchor_including_alignments(list_dict_alignments)
    table_anchor_paf_filtered = filter_table_with_anchor_alignment_position(table_anchor_paf, dict_region, anchor_min_align, anchor_min_align2hang, anchor_max_wronghang, anchor_min_idy)
    return table_anchor_paf_filtered

def make_table_of_anchor_including_alignments(list_dict_alignments):
    table = pd.DataFrame()

    for dict_values in list_dict_alignments:
        dict_row = deepcopy(dict_values["alignment"])
        dict_row["front_anchor"] = dict_values["anchor_including"]["front"]
        dict_row["tail_anchor"] = dict_values["anchor_including"]["tail"]
        dict_row["cigar"] = dict_values["cigar"]

        table = pd.concat([table, pd.DataFrame(pd.Series(dict_row)).T])
    
    for bool_col in ["front_anchor", "tail_anchor"]:
        table[bool_col] = table[bool_col].astype(bool)
    for int_col in ["qrylen", "qry_alignstart", "qry_alignend", "reflen", "ref_alignstart", "ref_alignend", "matchlen", "alignlen"]:
        table[int_col] = table[int_col].astype(int)
    table["qry_alignend"] = table["qry_alignend"] - 1
    table["ref_alignend"] = table["ref_alignend"] - 1

    return table

def filter_table_with_anchor_alignment_position(table_anchor, lowqv_info, min_align, min_align2hang, max_wronghang, min_idy):
    table_anchor_filt = table_anchor.copy()

    # Filter by alignment length of anchored region
    table_anchor_filt["front_anchor_alignlen"] = table_anchor_filt.apply(lambda row: calculate_alignlen_of_anchor(row["qry_alignstart"], row["qry_alignend"], lowqv_info["anchor_front"]["start"], lowqv_info["anchor_front"]["end"] if row["front_anchor"] else 0), axis = 1)
    table_anchor_filt["tail_anchor_alignlen"] = table_anchor_filt.apply(lambda row: calculate_alignlen_of_anchor(row["qry_alignstart"], row["qry_alignend"], lowqv_info["anchor_tail"]["start"], lowqv_info["anchor_tail"]["end"] if row["tail_anchor"] else 0), axis = 1)
    table_anchor_filt = table_anchor_filt[np.logical_or(~table_anchor_filt["front_anchor"], table_anchor_filt["front_anchor_alignlen"] >= min_align)]
    table_anchor_filt = table_anchor_filt[np.logical_or(~table_anchor_filt["tail_anchor"], table_anchor_filt["tail_anchor_alignlen"] >= min_align)]
    if table_anchor_filt.shape[0] == 0:
        return pd.DataFrame()
    
    # Filter by hanged alignment length of anchored region
    table_anchor_filt["front_hang"] = table_anchor_filt.apply(
        lambda row:
            sum(calculate_wronghang(row["qry_alignstart"], row["qry_alignend"], row["ref_alignstart"], row["ref_alignend"], row["reflen"], True, lowqv_info["anchor_front"]["start"], lowqv_info["anchor_front"]["end"], row["strand"])) if row["front_anchor"] else 0,
        axis = 1
    )
    table_anchor_filt["tail_hang"] = table_anchor_filt.apply(
        lambda row:
            sum(calculate_wronghang(row["qry_alignstart"], row["qry_alignend"], row["ref_alignstart"], row["ref_alignend"], row["reflen"], False, lowqv_info["anchor_tail"]["start"], lowqv_info["anchor_tail"]["end"], row["strand"])) if row["tail_anchor"] else 0,
        axis = 1
    )
    table_anchor_filt = table_anchor_filt[np.logical_and(table_anchor_filt["front_hang"] < max_wronghang, table_anchor_filt["tail_hang"] < max_wronghang)]
    if table_anchor_filt.shape[0] == 0:
        return pd.DataFrame()

    # Filter by align2hang ratio
    table_anchor_filt["front_align2hang"] = (table_anchor_filt["front_anchor_alignlen"] / (table_anchor_filt["front_anchor_alignlen"] + table_anchor_filt["front_hang"])).fillna(1)
    table_anchor_filt["tail_align2hang"] = (table_anchor_filt["tail_anchor_alignlen"] / (table_anchor_filt["tail_anchor_alignlen"] + table_anchor_filt["tail_hang"])).fillna(1)
    table_anchor_filt = table_anchor_filt[np.logical_and(table_anchor_filt["front_align2hang"] >= min_align2hang, table_anchor_filt["tail_align2hang"] >= min_align2hang)]
    if table_anchor_filt.shape[0] == 0:
        return pd.DataFrame()

    # Filter by Identity
    table_anchor_filt["front_identity"] = table_anchor_filt.apply(lambda row: calculate_identity_from_cigar_with_length(row["cigar"], row["front_anchor_alignlen"], True) if row["front_anchor"] else 1, axis = 1)
    table_anchor_filt["tail_identity"] = table_anchor_filt.apply(lambda row: calculate_identity_from_cigar_with_length(row["cigar"], row["tail_anchor_alignlen"], False) if row["tail_anchor"] else 1, axis = 1)
    table_anchor_filt = table_anchor_filt[np.logical_and(table_anchor_filt["front_identity"] >= min_idy, table_anchor_filt["tail_identity"] >= min_idy)]
    
    return table_anchor_filt

def calculate_alignlen_of_anchor(alignstart, alignend, anchor_start, anchor_end):
    set_positions_align = set(range(alignstart, alignend+1))
    set_positions_anchor = set(range(anchor_start, anchor_end+1))

    return len(set_positions_align & set_positions_anchor)

def calculate_wronghang(qry_alignstart, qry_alignend, ref_alignstart, ref_alignend, reflen, is_anchor_front, anchor_start, anchor_end, alignstrand):
    if is_anchor_front:
        if alignstrand == '+':
            lefthang = min(qry_alignstart - anchor_start, ref_alignstart)
            # lefthang = 0
            # righthang = 0 if qry_alignend >= anchor_end else anchor_end - qry_alignend
            righthang = 0 if qry_alignend >= anchor_end else np.inf
        else:
            # lefthang = 0 if qry_alignend >= anchor_end else anchor_end - qry_alignend
            lefthang = 0 if qry_alignend >= anchor_end else np.inf
            righthang = min(qry_alignstart - anchor_start, reflen - ref_alignend - 1)
            # righthang = 0
    else:
        if alignstrand == '+':
            # lefthang = 0 if qry_alignstart <= anchor_start else qry_alignstart - anchor_start
            lefthang = 0 if qry_alignstart <= anchor_start else np.inf
            righthang = min(anchor_end - qry_alignend, reflen - ref_alignend - 1)
            # righthang = 0
        else:
            lefthang = min(anchor_end - qry_alignend, ref_alignstart)
            # lefthang = 0
            # righthang = 0 if qry_alignstart <= anchor_start else qry_alignstart - anchor_start
            righthang = 0 if qry_alignstart <= anchor_start else np.inf
    return (lefthang, righthang)

def calculate_identity_from_cigar_with_length(cigar, length, is_front):
    list_cigar = elongate_cigar(cigar)
    if not is_front:
        list_cigar = list_cigar[::-1]

    list_cigar_of_desired = list()
    cnt = 0
    for cigar_char in list_cigar:
        list_cigar_of_desired.append(cigar_char)
        if cigar_char != 'D':
            cnt += 1
        if cnt == length:
            break
    
    count_cigar = Counter(list_cigar_of_desired)
    identity = count_cigar["="] / len(list_cigar_of_desired)
    return identity

def elongate_cigar(cigar):
    list_cigar_char = list()
    cigar_re = re.compile("([0-9]+)([A-Z,=]+)")
    for find_result in cigar_re.finditer(cigar):
        num_elong, cchar = find_result.groups()
        list_cigar_char.extend([cchar] * int(num_elong))
    return list_cigar_char

# %%
if __name__ == "__main__":
    import argparse

    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument("--path_paf")
    argument_parser.add_argument("--path_lowqv_region_info")
    argument_parser.add_argument("--path_output")
    argument_parser.add_argument("--anchor_min_idy", type = float, default=0.9)
    argument_parser.add_argument("--anchor_min_align", type = int, default=1000)
    argument_parser.add_argument("--anchor_max_wronghang", type = int, default=1000)
    argument_parser.add_argument("--anchor_min_align2hang", type = float, default=0.9)
    argument_parser.add_argument("--n_jobs", type = int, default=1)

    args = argument_parser.parse_args()

    run(
        args.path_paf,
        args.path_lowqv_region_info,
        args.path_output,
        args.anchor_min_idy,
        args.anchor_min_align,
        args.anchor_max_wronghang,
        args.anchor_min_align2hang,
        args.n_jobs
    )