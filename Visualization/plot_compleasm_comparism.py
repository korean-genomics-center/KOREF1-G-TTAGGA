#%%
import math

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec

#%%

dict_asm_to_compleasm = {
    "KOREF-S1-TTAGGA v1":{
        "Paternal":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Compleasm_Evaluation/Hap1/summary.txt",
        "Maternal":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Compleasm_Evaluation/Hap2/summary.txt" ,
        "Maternal w/o X":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Compleasm_Evaluation/Hap2.Autosome/summary.txt"   
    },
    "CHM13 v2":{
        "Haploid":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Compleasm_Evaluation/CHM13_v2/Haploid/summary.txt"
    },
    "T2T-YAO v1.1":{
        "Paternal":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Compleasm_Evaluation/T2T_YAO/Paternal/summary.txt",
        "Maternal":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Compleasm_Evaluation/T2T_YAO/Maternal/summary.txt"
    },
    "CN1 v1.0.1":{
        "Paternal":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Compleasm_Evaluation/CN1/Paternal/summary.txt",
        "Maternal":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Compleasm_Evaluation/CN1/Maternal/summary.txt"
    },
    "I002C v0.7":{
        "Paternal" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Compleasm_Evaluation/I002C/Paternal/summary.txt",
        "Maternal" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Compleasm_Evaluation/I002C/Maternal/summary.txt"
    }
}
# %%
def read_compleasm_summary(path_summary):
    dict_summary = dict()
    with open(path_summary, 'r') as fr:
        for line in fr:
            if line.startswith('#'):
                continue
            line = line.strip('\n')
            name = line.split(':')[0]
            count = int(line.split(':')[-1].split(',')[-1].strip())

            dict_summary[name] = count
    return dict_summary

dict_compleasm_result_for_asm = dict()

list_asmnametype = list()

dict_asmnametype_to_count = {
    'S':list(),
    'D':list(),
    'F':list(),
    'I':list(),
    'M':list(),
    'N':list()
}

for asmname in dict_asm_to_compleasm.keys():
    dict_compleasm_result_for_asm[asmname] = dict()
    for asmtype, path_summary in dict_asm_to_compleasm[asmname].items():
        dict_compleasm_result_for_asm[asmname][asmtype] = read_compleasm_summary(path_summary)

        list_asmnametype.append(f"{asmname} ({asmtype})")

        for key in dict_asmnametype_to_count.keys():
            dict_asmnametype_to_count[key].append(dict_compleasm_result_for_asm[asmname][asmtype][key])
# %%
dict_key_to_name = {
    'S':"Complete, Single-copy",
    'D':"Complete, Duplicated",
    'F':"Fragmented (Partial)",
    'I':"Fragmented (Chimeric)",
    'M':"Missing"
}

dict_key_to_color = {
    'S':"cornflowerblue",
    'D':"mediumblue",
    'F':"gold",
    'I':"darkorange",
    'M':"firebrick"
}

plt.rcParams["font.size"] = 14
plt.rcParams["font.family"] = "sans-serif"

plt.figure(figsize = (6, 4))

bottom = np.zeros(len(list_asmnametype))

for key in ['S', 'D', 'F', 'I', 'M']:
    values = dict_asmnametype_to_count[key]
    plt.barh(
        list_asmnametype[::-1], values[::-1],
        left = bottom,
        color = dict_key_to_color[key],
        label = dict_key_to_name[key]
    ),
    bottom += np.array(values[::-1])

# plt.xticks(rotation = -45, va = "center", ha = "left", rotation_mode="anchor")
plt.legend(loc = "center left", bbox_to_anchor = (1.01, 0.5))
# %%
