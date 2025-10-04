#%%
from copy import deepcopy

from joblib import Parallel, delayed
from zlib_ng import gzip_ng
import pandas as pd


def run(path_vcf, path_mosdepth_result, path_save, n_jobs = 2):
    table_vcf = read_vcf_as_table(path_vcf)
    table_vcf_sorted = table_vcf.sort_values(by = ["#CHROM", "POS"])
    dict_chrname_to_variant_pos = get_dict_of_variant_pos_on_table_vcf(table_vcf_sorted)
    
    dict_all_pos_to_depth = get_depth_of_vcf_pos(dict_chrname_to_variant_pos, path_mosdepth_result, n_jobs)
    
    table_depth = convert_dict_pos_to_depth_into_table(dict_chrname_to_variant_pos, dict_all_pos_to_depth,)
    
    table_depth.to_csv(path_save, sep = '\t', index = False)

def read_vcf_as_table(path_vcf):
    n_header = check_vcf_header_len(path_vcf)
    table_vcf = pd.read_csv(path_vcf, sep = '\t', skiprows = n_header)
    return table_vcf

def check_vcf_header_len(path_vcf):
    n_header = 0

    read_mode = 'rb' if path_vcf.endswith(".gz") else 'r'
    with open(path_vcf, read_mode) as fr:
        for line in fr:
            if hasattr(line, "decode"):
                line = line.decode()
            if line.startswith("#CHROM"):
                break
            else:
                n_header += 1
    return n_header

def get_depth_of_vcf_pos(dict_chrname_to_variant_pos, path_mosdepth_result, n_jobs):
    list_chrnames = dict_chrname_to_variant_pos.keys()
    with Parallel(n_jobs = min(n_jobs, len(dict_chrname_to_variant_pos))) as parallel:
        list_pos_to_depth = parallel(delayed(get_depth_of_vcf_pos_of_specific_chromosome)(
            dict_chrname_to_variant_pos[chrname],
            chrname,
            path_mosdepth_result
        ) for chrname in list_chrnames)
        
    dict_all_pos_to_depth = dict()
    for chrname, dict_pos_to_depth in zip(list_chrnames, list_pos_to_depth):
        dict_all_pos_to_depth[chrname] = dict_pos_to_depth
    return dict_all_pos_to_depth    

def get_depth_of_vcf_pos_of_specific_chromosome(list_pos_check, chrname_check, path_mosdepth_result):
    set_pos_check = set(list_pos_check)
    
    dict_pos_to_depth = dict()
    
    with gzip_ng.open(path_mosdepth_result, "rb") as frb:
        is_line_passed_chrname = False
        for line in frb:
            line = line.decode()
            list_values = line.strip('\n').split('\t')
            chrname = list_values[0]
            
            if chrname != chrname_check:
                if is_line_passed_chrname:
                    break
                continue
            else:
                is_line_passed_chrname = True
                
                range_start = int(list_values[1])+1
                range_end = int(list_values[2])
                
                for pos in list(set_pos_check):
                    if check_is_value_in_range((range_start, range_end), pos):
                        dict_pos_to_depth[pos] = int(list_values[3])
                        set_pos_check.remove(pos)
                if len(set_pos_check) == 0:
                    break
    return dict_pos_to_depth
    
def get_dict_of_variant_pos_on_table_vcf(table_vcf):
    dict_chrname_to_variant_pos = dict()
    for _, row in table_vcf.iterrows():
        if dict_chrname_to_variant_pos.get(row["#CHROM"]) == None:
            dict_chrname_to_variant_pos[row["#CHROM"]] = list()
        
        dict_chrname_to_variant_pos[row["#CHROM"]].append(row["POS"])
    return dict_chrname_to_variant_pos

def check_is_value_in_range(ranges, value):
    return value >= ranges[0] and value <= ranges[1]

def convert_dict_pos_to_depth_into_table(dict_chrname_to_variant_pos, dict_all_pos_to_depth):
    list_tables_concat = list()
    
    for chrname in dict_chrname_to_variant_pos.keys():
        table_depth_of_chrname = pd.DataFrame({"POS" : dict_chrname_to_variant_pos[chrname]})
        table_depth_of_chrname["#CHROM"] = chrname
        table_depth_of_chrname["Depth"] = table_depth_of_chrname["POS"].apply(lambda val: dict_all_pos_to_depth[chrname].get(val))
        
        list_tables_concat.append(table_depth_of_chrname[["#CHROM", "POS", "Depth"]].copy())
    return pd.concat(list_tables_concat)

#%%
if __name__ == "__main__":
    import argparse
    
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument("--vcf")
    argument_parser.add_argument("--mosdepth")
    argument_parser.add_argument("--out")
    argument_parser.add_argument("--n_jobs", type = int, default = 1)
    
    args = argument_parser.parse_args()
    
    run(
        args.vcf,
        args.mosdepth,
        args.out,
        args.n_jobs
    )