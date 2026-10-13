#%%
from pathlib import Path
import sys,gzip,math

from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

sys.path.append(str(Path(__file__).parents[2] / "common"))
from fasta_parser import read_fasta_as_dict


#%%
path_depth = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/GCI_Evaluation/GCI_Result/Hap2/GCI_Result.hap2.HiFi_ONT.Winnowmap2_Minimap2.triobin_m0a_two_type.depth"
path_issue = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/GCI_Evaluation/GCI_Result/Hap2/GCI_Result.hap2.HiFi_ONT.Winnowmap2_Minimap2.triobin_m0a_two_type.0.depth.bed"
path_fasta = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Final_Assembly/FinalASM.KOREF_S1.Trio_Specific.hap2.Reordered.fasta"

window_size = 10000

# %%
dict_fasta = read_fasta_as_dict(path_fasta)
# %%
list_chromosomes = list(dict_fasta.keys())
# %%
table_issue_bed = pd.read_csv(path_issue, sep = '\t', names = ["chr", "start", "end"])
# %%
dict_depth = dict()
with open(path_depth, "r") as fr:
    chrname = None
    
    for line in fr:
        line = line.strip('\n')
        if line.startswith('>'):
            chrname = line[1:]
            print(chrname)
            assert dict_depth.get(chrname) == None
            dict_depth[chrname] = list()
        else:
            val = int(line)
            dict_depth[chrname].append(val)
# %%
dict_depth_sum_window = dict()
dict_depth_n_window = dict()

for chrname in list_chromosomes:
    dict_depth_sum_window[chrname] = dict()
    dict_depth_n_window[chrname] = dict()
    
    list_depth = dict_depth[chrname]
    
    print(len(dict_fasta[chrname]), len(list_depth))
    
    for ind, depth in enumerate(list_depth):
        window_ind = math.floor(ind / window_size)
        
        if dict_depth_sum_window[chrname][window_ind] == None:
            dict_depth_sum_window[chrname][window_ind] = 0
            dict_depth_n_window[chrname][window_ind] = 0
        dict_depth_sum_window[chrname][window_ind] += depth
        dict_depth_n_window[chrname][window_ind] += 1

dict_depth_mean_window = dict()
for chrname in list_chromosomes:
    dict_depth_mean_window[chrname] = np.array(dict_depth_sum_window[chrname]) / np.array(dict_depth_n_window[chrname])
