#%%
import os
import pandas as pd
import numpy as np


# %%
def check_header_len(path_vcf):
    n_header = 0
    with open(path_vcf, 'r') as fr:
        for line in fr:
            if line.startswith("#CHROM"):
                break
            else:
                n_header += 1
    return n_header

def read_vcf_as_table(path_vcf):
    n_header = check_header_len(path_vcf)
    table_vcf = pd.read_csv(path_vcf, sep = '\t', skiprows = n_header)
    return table_vcf

# %%
def parse_info_data_from_vcf(info_data):
    list_data = info_data.split(';')
    
    list_flags = list()
    dict_item = dict()
    
    for data in list_data:
        if '=' in data:
            key, value = data.split('=')
            assert dict_item.get(key) == None
            dict_item[key] = value
        else:
            list_flags.append(data)
    return list_flags, dict_item


#%%
# MapQ20.Err0.05
# MapQ20.Err0.1

hap = 1

table_vcf = read_vcf_as_table(f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round1_SV_call/SV_Calling/Merged_SV/SV_Call.KOREF_GapclosedASM_Round1_SV.hap{hap}.hifi_ont.sametrio.filtered.iris.merged.split.homoalt_include.hetalt_only.filt_difftrio_homoref.filt_no_asmsv.vaf_0.5.minsupp_5.vcf")

# table_vcf = read_vcf_as_table(f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Align_Gapclosed_ASM/SV_Call.Inc_Amb/SV_Calling/Merged_SV/SV_Call.KOREF_GapclosedASM_SV_Polishing.Inc_Amb.hap{hap}.hifi_ont.sametrio.filtered.iris.merged.split.hetalt_only.filt_difftrio_homoref.filt_no_asmsv.vaf_0.5.minsupp_5.vcf")

series_info_flags = table_vcf["INFO"].apply(lambda val: parse_info_data_from_vcf(val)[0])
series_info_dict = table_vcf["INFO"].apply(lambda val: parse_info_data_from_vcf(val)[1])
# %%
table_vcf["SVtype"] = series_info_dict.apply(lambda dictval: dictval.get("SVTYPE"))
table_vcf["Length"] = series_info_dict.apply(lambda dictval: dictval.get("SVLEN"))
table_vcf["Coverage"] = series_info_dict.apply(lambda dictval: dictval.get("COVERAGE"))
table_vcf["Strand"] = series_info_dict.apply(lambda dictval: dictval.get("STRAND"))
table_vcf["Support"] = series_info_dict.apply(lambda dictval: dictval.get("SUPPORT"))
table_vcf["SUPP_VEC"] = series_info_dict.apply(lambda dictval: dictval.get("SUPP_VEC"))

#%%
table_vcf["Genotype_0"] = table_vcf["0_SAMPLE"].apply(lambda val: val.split(':')[0])
table_vcf["Genotype_1"] = table_vcf["1_SAMPLE"].apply(lambda val: val.split(':')[0])
table_vcf["DV_0"] = table_vcf["0_SAMPLE"].apply(lambda val: val.split(':')[3])
table_vcf["DV_1"] = table_vcf["1_SAMPLE"].apply(lambda val: val.split(':')[3])
table_vcf["DR_0"] = table_vcf["0_SAMPLE"].apply(lambda val: val.split(':')[4])
table_vcf["DR_1"] = table_vcf["1_SAMPLE"].apply(lambda val: val.split(':')[4])
#%%
table_vcf["VAF_0"] = table_vcf.apply(lambda row: None if row["DV_0"] == "NA" else int(row["DV_0"]) / (int(row["DV_0"]) + int(row["DR_0"])), axis = 1)
table_vcf["VAF_1"] = table_vcf.apply(lambda row: None if row["DV_1"] == "NA" else int(row["DV_1"]) / (int(row["DV_1"]) + int(row["DR_1"])), axis = 1)
#%%
table_vcf["Length_Calc"] = table_vcf.apply(
    lambda row: len(row["ALT"]) - len(row["REF"]) if set(row["ALT"]) <= {'A','T','G','C'} else None, axis = 1
)
#%%
is_ont = 0
#%%
from matplotlib import pyplot as plt

plt.rcParams['font.size'] = 20
plt.hist(table_vcf[table_vcf[f"DV_{is_ont}"] != "NA"][f"DV_{is_ont}"].astype(int), bins = 30, color = "gray")
plt.xlabel("Min Support")
plt.ylabel("# SV")
#%%
plt.rcParams['font.size'] = 20
plt.hist(table_vcf[table_vcf[f"DV_{is_ont}"] != "NA"][f"DV_{is_ont}"].astype(int), bins = np.linspace(0, 70, 15), color = "gray")
plt.xlabel("Min Support")
plt.ylabel("# SV")
#%%
plt.rcParams['font.size'] = 20
plt.hist(table_vcf[table_vcf[f"DV_{is_ont}"] != "NA"][f"DV_{is_ont}"].astype(int), bins = np.linspace(0, 20, 21), color = "gray")
plt.xlabel("Min Support")
plt.ylabel("# SV")
#%%
plt.hist(table_vcf[f"VAF_{is_ont}"].dropna(), bins = np.linspace(0, 1, 101), color = "gray")
plt.xlabel("VAF")
plt.ylabel("# SV")
plt.axvline(0.4)
#%%
from collections import Counter

print(Counter(table_vcf["SUPP_VEC"]))
# %%
print(Counter(list(zip(table_vcf["Genotype_0"], table_vcf["Genotype_1"]))))
# %%
#########################################
table_vcf_merfin = read_vcf_as_table(f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Align_Gapclosed_ASM/SV_Call.Inc_Amb/Merfin_Filtered/SV_Call.KOREF_GapclosedASM_SV_Polishing.Inc_Amb.hap{hap}.hifi_ont.sametrio.filtered.iris.merged.split.homoalt_include.hetalt_only.filt_difftrio_homoref.filt_no_asmsv.vaf_0.5.minsupp_5.call_merged.iris.polish.vcf")

list_merfin_svs = table_vcf_merfin.apply(lambda row: f"{row['#CHROM']}:{row['POS']}_{row['REF']}_{row['ALT']}", axis = 1)
table_vcf["SV_KEY"] = table_vcf.apply(lambda row: f"{row['#CHROM']}:{row['POS']}_{row['REF']}_{row['ALT']}", axis = 1)


table_vcf_merfin_in = table_vcf[table_vcf["SV_KEY"].isin(list_merfin_svs.to_list())]
# %%
import math
max_bin = math.ceil(table_vcf_merfin_in["Length_Calc"].abs().max()/500)*500
#%%
from matplotlib import pyplot as plt
import numpy as np

plt.rcParams['font.size'] = 20
plt.hist(table_vcf_merfin_in["Length_Calc"].abs(), bins = np.linspace(0, max_bin, int(max_bin/50+1)), color = "gray")
plt.xlabel("SV size (bp)")
plt.ylabel("# SVs")
plt.show()
# %%
from collections import Counter
Counter(table_vcf_merfin_in["SVtype"])
# %%
suppvec = "11"

from collections import Counter
from matplotlib import pyplot as plt
# %%
count_svtype_vcfin = Counter(table_vcf[table_vcf["SUPP_VEC"] == suppvec]["SVtype"])
count_svtype_vcfmerf = Counter(table_vcf_merfin_in[table_vcf_merfin_in["SUPP_VEC"] == suppvec]["SVtype"])

#%%
list_svtypes = sorted(count_svtype_vcfin.keys())

#%%
plt.figure(figsize = (1.5 * len(list_svtypes), 4))
plt.rcParams["font.size"] = 14
plt.bar(
    list(map(lambda val: val - 0.15, range(len(list_svtypes)))),
    list(map(lambda key: count_svtype_vcfin[key], list_svtypes)),
    color = "gray",
    width = 0.3,
    label = "Before Filtering"
)
for ind, key in enumerate(list_svtypes):
    plt.annotate(
        count_svtype_vcfin[key],
        (ind-0.15, count_svtype_vcfin[key]),
        ha = "center", va = "bottom"
    )
plt.bar(
    list(map(lambda val: val + 0.15, range(len(list_svtypes)))),
    list(map(lambda key: count_svtype_vcfmerf[key], list_svtypes)),
    color = "firebrick",
    width = 0.3,
    label = "Merfin Filtered"
)
for ind, key in enumerate(list_svtypes):
    plt.annotate(
        count_svtype_vcfmerf[key],
        (ind+0.15, count_svtype_vcfmerf[key]),
        ha = "center", va = "bottom"
    )

plt.xticks(list(range(len(list_svtypes))), list_svtypes)
plt.ylabel("# SVs")
plt.ylim(top = max(count_svtype_vcfin.values())*1.1)

for ind, key in enumerate(list_svtypes):
    plt.annotate(
        f"{count_svtype_vcfmerf[key]/count_svtype_vcfin[key]*100:.1f}%",
        (ind, max(count_svtype_vcfin.values())*1.11),
        ha = "center", va = "bottom", annotation_clip = False, fontsize = plt.rcParams["font.size"] + 2
    )
# plt.legend(loc = "center left", bbox_to_anchor = (1.01, 0.5))
plt.show()

# %%
from collections import Counter
from matplotlib import pyplot as plt
import numpy as np
list_chrnames = list(map(lambda val: f"chr{val}", range(1, 23))) + ["chrX", "chrY"]

count_vcfmerfin_chrname = Counter(table_vcf_merfin_in["#CHROM"])
count_vcfmerfin_ins_chrname = Counter(table_vcf_merfin_in[table_vcf_merfin_in["SVtype"] == "INS"]["#CHROM"])
count_vcfmerfin_del_chrname = Counter(table_vcf_merfin_in[table_vcf_merfin_in["SVtype"] == "DEL"]["#CHROM"])

plt.figure(figsize = (9, 3))
plt.rcParams['font.size'] = 20

arr_counts = np.zeros(len(list_chrnames))

plt.bar(list_chrnames, list(map(lambda key: count_vcfmerfin_ins_chrname[key], list_chrnames)), color = "cornflowerblue", bottom = arr_counts, label = "Insertion")

arr_counts += np.array(list(map(lambda key: count_vcfmerfin_ins_chrname[key], list_chrnames)))

plt.bar(list_chrnames, list(map(lambda key: count_vcfmerfin_del_chrname[key], list_chrnames)), color = "salmon", bottom = arr_counts, label = "Deletion")


plt.xlabel("")
plt.ylabel("# SVs")
plt.xticks(rotation = 90)
plt.xlim(-0.5, len(list_chrnames)-0.5)
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1], loc = "center left", bbox_to_anchor = (1.01, 0.5))

for ind, chrname in enumerate(list_chrnames):
    plt.annotate(count_vcfmerfin_chrname[chrname], (ind, count_vcfmerfin_chrname[chrname]), va = "bottom", ha = "center")
plt.ylim(top = max(count_vcfmerfin_chrname.values())*1.2)
plt.show()
# %%
