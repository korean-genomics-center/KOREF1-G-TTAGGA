#%%
import math

import pandas as pd
import numpy as np


def run(
    path_qv_per_window,
    path_save_prefix,
    qv_thresh = 40,
    extra_ranges_collect = 5
):
    path_save_bed = f"{path_save_prefix}.bed"
    path_save_front_anchor_bed = f"{path_save_prefix}.front_anchor.bed"
    path_save_tail_anchor_bed = f"{path_save_prefix}.tail_anchor.bed"
    path_save_info = f"{path_save_prefix}.tsv"

    table_qv = pd.read_csv(path_qv_per_window, names = ["seqname", "wrong", "total", "qv", "error"], sep = '\t')

    table_qv.loc[np.isinf(table_qv["qv"]), "qv"] = 99
    table_qv["chrname"] = table_qv["seqname"].apply(lambda val: val.split(':')[-2])
    table_qv["range_start"] = table_qv["seqname"].apply(lambda val: int(val.split(':')[-1].split('-')[0]))
    table_qv["range_end"] = table_qv["seqname"].apply(lambda val: int(val.split(':')[-1].split('-')[1]))
    table_qv["range_mid"] = (table_qv["range_start"]+table_qv["range_end"]) / 2

    list_regions_get = list()
    list_regions_front_anchor = list()
    list_regions_tail_anchor = list()

    file_write_info = open(path_save_info, 'w')
    file_write_info.write(get_context_to_write(["chr", "range_start", "range_end", "lowqv_start", "lowqv_end"]))

    for chrname in table_qv["chrname"].unique():
        table_qv_chr = table_qv[table_qv["chrname"] == chrname].reset_index(drop = True)

        chr_start = 0
        chr_end = table_qv_chr.shape[0]-1

        table_qv_chr_low = table_qv_chr[table_qv_chr["qv"] < qv_thresh]
        if table_qv_chr_low.shape[0] == 0:
            continue
        
        list_ind_qv_low = list(table_qv_chr_low.index)
        list_inds_ranged = list()

        for ind_low in list_ind_qv_low:
            if ind_low == chr_start or ind_low == chr_end:
                continue
            if ind_low in list_inds_ranged:
                continue

            list_inds_for_range = [ind_low]

            ind_low_start = ind_low
            ind_low_end = ind_low

            # extend front
            n_extend_front = 0
            ind_add_front = ind_low-1
            while 1:
                if ind_add_front == 0:
                    list_inds_for_range.append(ind_add_front)
                    break
                if table_qv_chr.loc[ind_add_front, "qv"] > 80:
                    list_inds_for_range.append(ind_add_front)
                    n_extend_front += 1
                else:
                    ind_low_start = ind_add_front
                    list_inds_for_range.append(ind_add_front)
                    n_extend_front = 0
                ind_add_front -= 1
                if n_extend_front == extra_ranges_collect:
                    break
            
            # extend tail
            n_extend_tail = 0
            ind_add_tail = ind_low+1
            while 1:
                if ind_add_tail == chr_end:
                    list_inds_for_range.append(ind_add_tail)
                    break
                if table_qv_chr.loc[ind_add_tail, "qv"] > 80:
                    list_inds_for_range.append(ind_add_tail)
                    n_extend_tail += 1
                else:
                    ind_low_end = ind_add_tail
                    list_inds_for_range.append(ind_add_tail)
                    n_extend_tail = 0
                ind_add_tail += 1
                if n_extend_tail == extra_ranges_collect:
                    break
            
            list_inds_for_range = sorted(list_inds_for_range)
            list_inds_ranged.extend(list_inds_for_range)
            list_regions_get.append(
                (
                    chrname, table_qv_chr.loc[list_inds_for_range[0], "range_start"], table_qv_chr.loc[list_inds_for_range[-1], "range_end"]
                )
            )
            list_regions_front_anchor.append(
                (
                    chrname, table_qv_chr.loc[list_inds_for_range[0], "range_start"], table_qv_chr.loc[list_inds_for_range[extra_ranges_collect-1], "range_end"]
                )
            )
            list_regions_tail_anchor.append(
                (
                    chrname, table_qv_chr.loc[list_inds_for_range[-extra_ranges_collect], "range_start"], table_qv_chr.loc[list_inds_for_range[-1], "range_end"]
                )
            )

            file_write_info.write(
                get_context_to_write(
                    [
                        chrname, 
                        table_qv_chr.loc[list_inds_for_range[0], "range_start"], 
                        table_qv_chr.loc[list_inds_for_range[-1], "range_end"], 
                        table_qv_chr.loc[ind_low_start, "range_start"], 
                        table_qv_chr.loc[ind_low_end, "range_end"]
                    ]
                )
            )
    file_write_info.close()

    with open(path_save_bed, 'w') as fw:
        for region in list_regions_get:
            fw.write('\t'.join(list(map(str, region))) + '\n')
    
    with open(path_save_front_anchor_bed, 'w') as fw:
        for region in list_regions_front_anchor:
            fw.write('\t'.join(list(map(str, region))) + '\n')

    with open(path_save_tail_anchor_bed, 'w') as fw:
        for region in list_regions_tail_anchor:
            fw.write('\t'.join(list(map(str, region))) + '\n')


def get_context_to_write(list_values):
    return '\t'.join(list(map(str, list_values))) + '\n'
#%%
if __name__ == "__main__":
    import argparse

    argument_parser = argparse.ArgumentParser()

    argument_parser.add_argument("--path_qv")
    argument_parser.add_argument("--out_prefix")
    argument_parser.add_argument("--qv_thresh", type = int)
    argument_parser.add_argument("--extra_ranges", type = int)
    
    args = argument_parser.parse_args()

    run(
        args.path_qv,
        args.out_prefix,
        args.qv_thresh,
        args.extra_ranges
    )