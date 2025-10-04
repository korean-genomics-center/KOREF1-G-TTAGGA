#%%
import os
import pandas as pd


# %%
def check_header_len(path_vcf):
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

table_vcf = read_vcf_as_table(f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round4_SV_Call/SV_Call_Merged/Merged_Call/SV_Call.KOREF_PolishedASM_Round4_SV.hap2.hifi_ont.sametrio.filtered.iris.final_evaluated.indel_only.additional_support_of_unclassified_reads_on_depth_5.vcf")

col_sample = "0_SAMPLE__1_SAMPLE"

series_info_flags = table_vcf["INFO"].apply(lambda val: parse_info_data_from_vcf(val)[0])
series_info_dict = table_vcf["INFO"].apply(lambda val: parse_info_data_from_vcf(val)[1])
#%%
table_vcf["SVLEN"] = table_vcf.apply(lambda row: len(row["ALT"]) - len(row["REF"]), axis = 1)

table_vcf_manual_check = table_vcf[table_vcf["SVLEN"].abs() > 500].sort_values(by = ["#CHROM", "POS"])[["#CHROM", "POS", "ID", col_sample, "SVLEN"]]


for val in list(table_vcf_manual_check["ID"]):
    print(f""" "{val}" :  """)

table_vcf_manual_check
# %%
table_vcf["SVtype"] = series_info_dict.apply(lambda dictval: dictval.get("SVTYPE"))
table_vcf["Length"] = series_info_dict.apply(lambda dictval: dictval.get("SVLEN"))
table_vcf["Coverage"] = series_info_dict.apply(lambda dictval: dictval.get("COVERAGE"))
table_vcf["Strand"] = series_info_dict.apply(lambda dictval: dictval.get("STRAND"))
table_vcf["Genotype"] = table_vcf[col_sample].apply(lambda val: val.split(':')[0])
#%%
table_vcf["Support"] = series_info_dict.apply(lambda dictval: int(dictval.get("SUPPORT")))
#%%
table_vcf["Length_Calc"] = table_vcf.apply(
    lambda row: len(row["ALT"]) - len(row["REF"]) if set(row["ALT"]) <= {'A','T','G','C'} else None, axis = 1
)

# %%
table_vcf["Ref_Read"] = table_vcf.apply(lambda row: int(dict(zip(row["FORMAT"].split(':'), row[col_sample].split(':'))).get("DR")), axis = 1)
table_vcf["Var_Read"] = table_vcf.apply(lambda row: int(dict(zip(row["FORMAT"].split(':'), row[col_sample].split(':'))).get("DV")), axis = 1)
# %%
table_vcf["Var_Ratio"] = table_vcf["Var_Read"] / (table_vcf["Var_Read"]+table_vcf["Ref_Read"])
# %%
from collections import Counter
from matplotlib import pyplot as plt
plt.rcParams["font.size"] = 15
#%%
count_vartype = Counter(table_vcf["SVtype"])
count_homhet = Counter(table_vcf["Genotype"])

#%%
def plot_counter(count_data, keys_sorted = None):
    if keys_sorted == None:
        keys_sorted = sorted(count_data.keys(), key = count_data.__getitem__, reverse = True)

    plt.bar(
        keys_sorted,
        list(map(count_data.__getitem__, keys_sorted)),
        color = "gray"
    )
    plt.show()

# %%
plt.figure(figsize = (1.5, 4))
plot_counter(count_vartype)
# %%
plt.figure(figsize = (1.5, 4))
plot_counter(count_homhet, ["0/1", "1/1"])

# %%
import numpy as np

plt.figure(figsize = (5, 4))
plt.hist(table_vcf[table_vcf["Genotype"] == "1/1"]["Var_Ratio"], bins = np.linspace(0, 1, 101), color = "gray")
plt.xlabel("Variant Allele Frequency")
plt.ylabel("# SV")
# %%
plt.figure(figsize = (5, 4))
plt.hist(table_vcf[table_vcf["Genotype"] == "0/1"]["Var_Ratio"], bins = np.linspace(0, 1, 101), color = "gray")
plt.xlabel("Variant Allele Frequency")
plt.ylabel("# SV")
#%%
import math

plt.figure(figsize = (5, 4))

max_supp = table_vcf["Support"].max()
xlim_right = math.ceil(max_supp / 20) * 20

plt.hist(table_vcf["Support"], color = "gray", bins = np.linspace(0, xlim_right, xlim_right+1))
plt.xlabel("# Support")
plt.ylabel("# SVs")
plt.xlim(right  = 100)
# %%
import math
import numpy as np

max_svsize = table_vcf["Length_Calc"].dropna().abs().max()
bin_max = math.ceil(max_svsize / 1000) * 1000

plt.figure(figsize = (5, 4))
plt.hist(table_vcf["Length_Calc"].dropna().abs(), color = "gray", bins = np.linspace(0, bin_max, int(bin_max/400)+1))
plt.xlabel("Size of SV (bp)")
plt.ylabel("# SV")
# plt.xlim(0, bin_max)
plt.show()
#%%

plt.figure(figsize = (5, 4))
plt.hist(table_vcf[table_vcf["Length_Calc"] < 0]["Length_Calc"].dropna().abs(), color = "gray", bins = np.linspace(0, bin_max, int(bin_max/200)+1))
plt.xlabel("Size of Deletion SV (bp)")
plt.ylabel("# Deletion SV")
# plt.xlim(0, bin_max)
# plt.yscale("log")
plt.show()

#%%

plt.figure(figsize = (5, 4))
plt.hist(table_vcf[table_vcf["Length_Calc"] > 0]["Length_Calc"].dropna().abs(), color = "gray", bins = np.linspace(0, bin_max, int(bin_max/200)+1))
plt.xlabel("Size of Insertion SV (bp)")
plt.ylabel("# Insertion SV")
# plt.xlim(0, bin_max)
# plt.yscale("log")
plt.show()
#%%
import numpy as np

print("Mean", np.mean(table_vcf["Length_Calc"].dropna().abs()))
print("Max", np.max(table_vcf["Length_Calc"].dropna().abs()))
print("Median", np.median(table_vcf["Length_Calc"].dropna().abs()))
# %%
table_vcf_notna = table_vcf[table_vcf["Length_Calc"].notna()]

table_vcf_notna.sort_values(by = ["Length_Calc"])
# %%
table_vcf_notna["SVtype"] = table_vcf_notna["Length_Calc"].apply(lambda val: "INS" if val > 0 else "DEL")
# %%
from collections import Counter
from matplotlib import pyplot as plt
import numpy as np
list_chrnames = list(map(lambda val: f"chr{val}", range(1, 23))) + ["chrX", "chrY"]

count_vcfmerfin_chrname = Counter(table_vcf_notna["#CHROM"])
count_vcfmerfin_ins_chrname = Counter(table_vcf_notna[table_vcf_notna["SVtype"] == "INS"]["#CHROM"])
count_vcfmerfin_del_chrname = Counter(table_vcf_notna[table_vcf_notna["SVtype"] == "DEL"]["#CHROM"])

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