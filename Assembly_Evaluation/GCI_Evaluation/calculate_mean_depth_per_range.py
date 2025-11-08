#%%
import numpy as np

import os, gzip


def main(path_depth, path_save, window_size = 50000):
    dict_meandepth = dict()
    
    is_gzip = False
    if path_depth.endswith(".gz"):
        file_reader = gzip.open(path_depth, "rb")
        is_gzip = True
            
    else:
        file_reader = open(path_depth, "r")
    
    chrname = None
    list_depth_tmp = list()
    
    for line in file_reader:
        if is_gzip:
            line = line.decode()
        line = line.strip('\n')     
        if line.startswith('>'):
            if len(list_depth_tmp) != 0:
                dict_meandepth[chrname].append(np.mean(list_depth_tmp))
            
            chrname = line[1:]
            assert dict_meandepth.get(chrname) == None, chrname
            
            dict_meandepth[chrname] = list()
            list_depth_tmp = list()
        else:
            depth = int(line)
            list_depth_tmp.append(depth)
            if len(list_depth_tmp) == window_size:
                dict_meandepth[chrname].append(np.mean(list_depth_tmp))
                list_depth_tmp = list()
    if len(list_depth_tmp) != 0:
        dict_meandepth[chrname].append(np.mean(list_depth_tmp))
    
    
    with open(path_save, 'w') as fw:
        for chrname, list_depth in dict_meandepth.items():
            fw.write(f">{chrname}\n")
            for depth in list_depth:
                fw.write(f"{depth}\n")
    
if __name__ == "__main__":
    from joblib import Parallel, delayed
    
    window_size = 50000
    
    dict_hap_to_parental = {
        1:'p',
        2:'m'
    }
    
    list_combinations = list()
    for hap in [1, 2]:
        for data in ["hifi", "nano", "two_type"]:
            list_combinations.append((hap, data))
    
    with Parallel(n_jobs = len(list_combinations)) as parallel:
        parallel(delayed(main)(
            f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/GCI_Evaluation/GCI_Result/Hap{hap}/GCI_Result.hap{hap}.HiFi_ONT.Winnowmap2_Minimap2.triobin_{dict_hap_to_parental[hap]}0a_{data}.depth.gz",
            f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/GCI_Evaluation/GCI_Result/Hap{hap}/GCI_Result.hap{hap}.HiFi_ONT.Winnowmap2_Minimap2.triobin_{dict_hap_to_parental[hap]}0a_{data}.depth.mean_{window_size}.gz",
            window_size = window_size
        ) for hap, data in list_combinations)