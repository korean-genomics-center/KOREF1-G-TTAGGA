#%%
import os
import pandas as pd
from Bio.bgzf import BgzfReader

# %%
def check_header_len(path_vcf):
    n_header = 0

    if path_vcf.endswith(".gz"):
        with BgzfReader(path_vcf, 'rb') as fr:
            for line in fr:
                line = line.decode()
                if line.startswith("#CHROM"):
                    break
                else:
                    n_header += 1
    else:
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

hap = 2
# triobin = 'pm0a'
# tgs = "HiFi"

# path_vcf = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Align_Polished_ASM.From_R1.SV/SNV_Call/SNV_Calling.{tgs}.Hap{hap}.Triobin_{triobin}/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.R1_Sv.hap{hap}.{tgs}.triobin_{triobin}.PASS.NonHomRef.biallelic.GQ30.VAF05.Norm.vcf.gz"

path_vcf = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Align_Polished_ASM.From_R1.SV/Polishing_Target/Merfin_Filtered/HiFi.Hap{hap}.Triobin_pm0a/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.R1_Sv.hap{hap}.HiFi.triobin_pm0a.PASS.NonHomRef.both_biallelic_multiallelic.biallelic_VAF04.ForPolishing.sorted.polish.vcf"
#%%
table_vcf = read_vcf_as_table(path_vcf)

series_info_flags = table_vcf["INFO"].apply(lambda val: parse_info_data_from_vcf(val)[0])
series_info_dict = table_vcf["INFO"].apply(lambda val: parse_info_data_from_vcf(val)[1])
#%%
series_sample_dict = table_vcf[["FORMAT", "Sample"]].apply(lambda row: dict(zip(row["FORMAT"].split(':'), row["Sample"].split(':'))), axis = 1)
# %%
table_vcf["GT"] = list(map(lambda dict_val: dict_val["GT"], series_sample_dict))
table_vcf["VAF"] = list(map(lambda dict_val: dict_val["VAF"], series_sample_dict))
table_vcf["GQ"] = list(map(lambda dict_val: dict_val["GQ"], series_sample_dict))
table_vcf["DP"] = list(map(lambda dict_val: dict_val["DP"], series_sample_dict))
#%%
table_vcf["Is_Biallelic"] = table_vcf["ALT"].apply(lambda val: len(val.split(',')) == 1)
#%%
table_vcf["GQ"] = table_vcf["GQ"].astype(int)
table_vcf["DP"] = table_vcf["DP"].astype(int)

table_vcf["N_ALT"] = table_vcf["ALT"].apply(lambda val: len(val.split(',')))

#%%
from matplotlib import pyplot as plt
import seaborn as sns
import math
import numpy as np
plt.rcParams["font.size"] = 15
list_n_alts = sorted(table_vcf["N_ALT"].unique())

fig, ax = plt.subplots(len(list_n_alts), 1, figsize = (7, 2*len(list_n_alts)))
axes = ax.flatten()

colors = sns.color_palette("Set2", len(list_n_alts))
max_gq = int(math.ceil(table_vcf["GQ"].max()/10)*10* 1.2) 
all_handles = list()
all_labels = list()
for ind, n_alt in enumerate(list_n_alts):
    histogram = axes[ind].hist(table_vcf[table_vcf["N_ALT"] == n_alt]["GQ"], bins = np.linspace(1, max_gq, max_gq+1), color = colors[ind], label = f"# ALT: {n_alt}")
    axes[ind].axvline(table_vcf[table_vcf["N_ALT"] == n_alt]["GQ"].max()+1, linewidth = 1.5, linestyle = "--", color = "firebrick")
    axes[ind].annotate(f"Max GQ:\n{table_vcf[table_vcf['N_ALT'] == n_alt]['GQ'].max()}", (table_vcf[table_vcf["N_ALT"] == n_alt]["GQ"].max()+2, 0), va = "bottom", ha = "left")
    handles_axes, labels_axes = axes[ind].get_legend_handles_labels()
    all_handles.append(handles_axes[0])
    all_labels.append(labels_axes[0])
axes[ind].set_xlabel("Genotype Quality")
plt.legend(all_handles, all_labels, loc = "center left", bbox_to_anchor = (1.01, 0.5))
plt.show()
# %%
from collections import Counter
from matplotlib import pyplot as plt
plt.rcParams["font.size"] = 15
#%%
count_is_biallelic = Counter(table_vcf["Is_Biallelic"])
#%%
def plot_counter(count_data, keys_sorted = None):
    if keys_sorted == None:
        keys_sorted = sorted(count_data.keys(), key = count_data.__getitem__, reverse = True)

    plt.bar(
        list(map(str, keys_sorted)),
        list(map(count_data.__getitem__, keys_sorted)),
        color = "gray"
    )

# %%
def convert_number_to_kmg(val):
    if val >= 1e9:
        val_annot = f"{val/1e9:.1f}G"
    elif val >= 1e6:
        val_annot = f"{val/1e6:.1f}M"
    elif val >= 1e3:
        val_annot = f"{val/1e3:.1f}K"
    else:
        val_annot = f"{val}"
    return val_annot
#%%
plt.figure(figsize = (1.5, 4))
plot_counter(count_is_biallelic, keys_sorted = [True, False])
plt.xlabel("Is Biallelic")
for ind, key in enumerate([True, False]):
    plt.annotate(convert_number_to_kmg(count_is_biallelic[key]), (ind, count_is_biallelic[key]), va = "bottom", ha = "center")
plt.ylim(top = max(count_is_biallelic.values())*1.1)
# plt.xlim(-0.75, 1.75)
plt.show()
# %%

# %%
import numpy as np

plt.figure(figsize = (5, 4))
plt.hist(table_vcf["GQ"], bins = np.linspace(0, 100, 101), color = "gray")
plt.xlabel("Genotype Quality")
plt.ylabel("# SNV")
# %%
import numpy as np
import math

dp_max = table_vcf["DP"].max()
dp_max_ceil = math.ceil(dp_max / 100) * 100

plt.figure(figsize = (5, 4))
plt.hist(table_vcf["DP"], bins = np.linspace(0, dp_max_ceil, 101), color = "gray")
plt.xlabel("Depth")
plt.ylabel("# SNV")
plt.show()
# %%
import numpy as np
import math

plt.figure(figsize = (5, 4))
plt.hist(table_vcf["DP"], bins = np.linspace(0, 100, 101), color = "gray")
plt.xlabel("Depth")
plt.ylabel("# SNV")
plt.show()
# %%
import numpy as np

plt.figure(figsize = (5, 4))
plt.hist(table_vcf[table_vcf["Is_Biallelic"]]["GQ"], bins = np.linspace(0, 100, 101), color = "gray")
plt.xlabel("Genotype Quality")
plt.ylabel("# SNV")
plt.show()
# %%
import numpy as np

plt.figure(figsize = (5, 4))
plt.hist(table_vcf[~table_vcf["Is_Biallelic"]]["GQ"], bins = np.linspace(0, 100, 101), color = "gray")
plt.xlabel("Genotype Quality")
plt.ylabel("# SNV")
plt.show()

# %%

#%%
#################################################
table_vcf_bi = table_vcf[table_vcf["Is_Biallelic"]]
table_vcf_bi["VAF"] = table_vcf_bi["VAF"].astype(float)

#%%
table_vcf_bi["Length_Calc"] = table_vcf_bi.apply(
    lambda row: len(row["ALT"]) - len(row["REF"]) if set(row["ALT"]) <= {'A','T','G','C'} else None, axis = 1
)

# %%
def get_snv_type_from_len(length):
    if length == 0:
        return "SNV"
    elif length < 0:
        return "DEL"
    else:
        return "INS"

table_vcf_bi["SNV_Type"] = table_vcf_bi["Length_Calc"].apply(get_snv_type_from_len)
# %%
count_bi_gt = Counter(table_vcf_bi["GT"])
count_bi_type = Counter(table_vcf_bi["SNV_Type"])
#%%
plt.figure(figsize = (1.5, 4))
plot_counter(count_bi_gt, keys_sorted = ["0/1", "1/1"])
for ind, key in enumerate(["0/1", "1/1"]):
    plt.annotate(convert_number_to_kmg(count_bi_gt[key]), (ind, count_bi_gt[key]), va = "bottom", ha = "center")
plt.ylim(top = max(count_bi_gt.values())*1.1)
# plt.xlim(-0.75, 1.75)
plt.show()
#%%
plt.figure(figsize = (3.5, 4))
plot_counter(count_bi_type, keys_sorted = ["SNV", "INS", "DEL"])
for ind, key in enumerate(["SNV", "INS", "DEL"]):
    plt.annotate(convert_number_to_kmg(count_bi_type[key]), (ind, count_bi_type[key]), va = "bottom", ha = "center")
plt.ylim(top = max(count_bi_type.values())*1.1)
plt.xlim(-0.75, 2.75)
plt.show()
# %%
import numpy as np

plt.figure(figsize = (5, 4))
plt.hist(table_vcf_bi[table_vcf_bi["GT"] == "0/1"]["VAF"], bins = np.linspace(0, 1, 101), color = "gray")
plt.xlabel("Variant Allele Frequency (VAF)")
plt.ylabel("# SNV")
# %%
import numpy as np

plt.figure(figsize = (5, 4))
plt.hist(table_vcf_bi[table_vcf_bi["GT"] == "1/1"]["VAF"], bins = np.linspace(0, 1, 101), color = "gray")
plt.xlabel("Variant Allele Frequency (VAF)")
plt.ylabel("# SNV")
# %%
# %%
import numpy as np

plt.figure(figsize = (5, 4))
plt.hist(table_vcf_bi["GQ"], bins = np.linspace(0, 100, 101), color = "gray")
plt.xlabel("Genotype Quality (GQ)")
plt.ylabel("# SNV")

# %%
import numpy as np

# table_draw = table_vcf_bi
table_draw = table_vcf_bi[table_vcf_bi["GT"] == "0/1"]
# table_draw = table_vcf_bi[table_vcf_bi["GT"] == "1/1"]
# table_draw = table_draw[table_draw["GQ"] > 30]

plt.figure(figsize = (5, 4))
plt.scatter(table_draw["VAF"], table_draw["GQ"], s = 3, color = "gray")
plt.axhline(10, linewidth = 1.5, color = "firebrick", linestyle = "--")
plt.xlabel("Variant Allele Frequency")
plt.ylabel("Genotype Quality")
plt.show()
# %%# %%
import numpy as np

table_draw = table_vcf_bi[table_vcf_bi["GQ"] > 10]

plt.figure(figsize = (5, 4))
plt.hist(table_draw[table_draw["GT"] == "1/1"]["VAF"], bins = np.linspace(0, 1, 101), color = "gray")
plt.xlabel("Variant Allele Frequency (VAF)")
plt.ylabel("# SNV")

# %%

#%%
#################################################
table_vcf_mul = table_vcf[~table_vcf["Is_Biallelic"]]
# %%
list_gt_mul = list(table_vcf_mul["GT"].unique())
list_gt_mul_het = list(filter(lambda val: '0' in val, list_gt_mul))
list_gt_mul_altonly = list(set(list_gt_mul) - set(list_gt_mul_het))
# %%
table_vcf_mul["Is_NoRef"] = table_vcf_mul["GT"].isin(list_gt_mul_altonly)
#%%
table_vcf_mul["Is_Het"] = table_vcf_mul["GT"].apply(lambda val: len(set(val.split('/'))) > 1)
#%%
table_vcf_mul["Max_VAF"] = table_vcf_mul["VAF"].apply(lambda val: max(list(map(float, val.split(',')))))
#%%
plt.figure(figsize = (5, 4))
plt.hist(table_vcf_mul["Max_VAF"], bins = np.linspace(0, 1, 101), color = "gray")
plt.xlabel("Max VAF of Multiallelic Variants")
plt.ylabel("# SNV")
plt.axvline(0.5, linewidth = 1.5, linestyle = "--", color = "firebrick")
plt.show()
# %%
table_vcf_mul_nonrefhet = table_vcf_mul[np.logical_and(table_vcf_mul["Is_NoRef"], table_vcf_mul["Is_Het"])]
# %%
table_vcf_mul_nonrefhet["VAF_1"] = table_vcf_mul_nonrefhet.apply(
    lambda row: max(
        list(map(
            lambda ind: list(map(float, row["VAF"].split(',')))[ind],
            list(map(lambda val: int(val)-1, row["GT"].split('/')))
        ))
    ),
    axis = 1
)
table_vcf_mul_nonrefhet["VAF_2"] = table_vcf_mul_nonrefhet.apply(
    lambda row: min(
        list(map(
            lambda ind: list(map(float, row["VAF"].split(',')))[ind],
            list(map(lambda val: int(val)-1, row["GT"].split('/')))
        ))
    ),
    axis = 1
)
# %%
import numpy as np

plt.figure(figsize = (5, 4))
plt.hist(table_vcf_mul[~table_vcf_mul["Is_NoRef"]]["GQ"], bins = np.linspace(0, 40, 41), color = "gray")
plt.xlabel("Genotype Quality")
plt.ylabel("# SNV")
plt.show()
#%%
import numpy as np

plt.figure(figsize = (5, 4))
plt.hist(table_vcf_mul[np.logical_and(table_vcf_mul["Is_NoRef"], ~table_vcf_mul["Is_Het"])]["GQ"], bins = np.linspace(0, 40, 41), color = "gray")
plt.xlabel("Genotype Quality")
plt.ylabel("# SNV")
plt.show()
#%%
import numpy as np

plt.figure(figsize = (5, 4))
plt.hist(table_vcf_mul[np.logical_and(table_vcf_mul["Is_NoRef"], table_vcf_mul["Is_Het"])]["GQ"], bins = np.linspace(0, 40, 41), color = "gray")
plt.axvline(25, linewidth = 1.5, linestyle = "--", color = "firebrick")
plt.xlabel("Genotype Quality")
plt.ylabel("# SNV")
plt.show()

# %%
plt.figure(figsize = (5, 4))
plt.scatter(table_vcf_mul_nonrefhet["VAF_2"], table_vcf_mul_nonrefhet["GQ"], color = "gray", s = 3)
plt.axhline(25, linewidth = 1.5, linestyle = "--", color = "firebrick")
plt.axvline(0.3, linewidth = 1.5, linestyle = "--", color = "firebrick")
plt.xlabel("Second Variant Allele Frequency")
plt.ylabel("Genotype Quality")
plt.show()
#%%
table_vcf_mul_nonrefhet_draw = table_vcf_mul_nonrefhet[table_vcf_mul_nonrefhet["GQ"] > 25]

plt.figure(figsize = (5, 4))
# plt.hist(table_vcf_mul_nonrefhet_draw["VAF_12"], bins = np.linspace(0, 1, 101), color = "gray", alpha = 0.5)
plt.hist(table_vcf_mul_nonrefhet_draw["VAF_1"], bins = np.linspace(0, 1, 101), color = "r", alpha = 0.5, label = "VAF A")
plt.hist(table_vcf_mul_nonrefhet_draw["VAF_2"], bins = np.linspace(0, 1, 101), color = "b", alpha = 0.5, label = "VAF B")
plt.axvline(0.3, linewidth = 1.5, linestyle = "--", color = "firebrick")
plt.xlabel("ALT Allele Frequency")
plt.ylabel("# SNV")
plt.legend(loc = "upper right", bbox_to_anchor = (0.95, 0.95))
plt.show()
#%%
plt.figure(figsize = (5, 4))
plt.hist(table_vcf_mul_nonrefhet[table_vcf_mul_nonrefhet["GQ"] > 25]["VAF_2"], bins = np.linspace(0, 1, 101), color = "gray")
plt.axvline(0.3, linewidth = 1.5, linestyle = "--", color = "firebrick")
plt.xlabel("Second Highest ALT Allele Frequency")
plt.ylabel("# SNV")
plt.show()
# %%
table_vcf_mul_nonrefhet["VAF_12"] = table_vcf_mul_nonrefhet["VAF_1"] + table_vcf_mul_nonrefhet["VAF_2"]
# %%
plt.figure(figsize = (5, 4))
plt.hist(table_vcf_mul_nonrefhet[table_vcf_mul_nonrefhet["GQ"] > 25]["VAF_12"], bins = np.linspace(0, 1, 101), color = "gray")
plt.axvline(0.7, linewidth = 1.5, linestyle = "--", color = "firebrick")
plt.xlabel("Sum of Top2 Alt Allele Frequency")
plt.ylabel("# SNV")
plt.show()
#%%
import seaborn as sns
plt.figure(figsize = (5, 4))
# sns.scatterplot(
#     data = table_vcf_mul_nonrefhet[table_vcf_mul_nonrefhet["GQ"] > 25],
#     x = "VAF_2",
#     y = "VAF_12",
#     hue = "GQ",
#     s = 5
# )
plt.scatter(
    table_vcf_mul_nonrefhet[table_vcf_mul_nonrefhet["GQ"] <= 25]["VAF_2"],
    table_vcf_mul_nonrefhet[table_vcf_mul_nonrefhet["GQ"] <= 25]["VAF_12"],
    s = 3,
    color = "gray"
)

plt.xlabel("VAF(B)")
plt.ylabel("VAF(A) + VAF(B)")
plt.axhline(0.8, linewidth = 1.5, color = "firebrick", linestyle = "--")
plt.axvline(0.3, linewidth = 1.5, color = "firebrick", linestyle = "--")
plt.show()
# %%
table_vcf_mul_nonrefhet_final = table_vcf_mul_nonrefhet[
    np.logical_and(
        np.logical_and(
            table_vcf_mul_nonrefhet["VAF_12"] > 0.7,
            table_vcf_mul_nonrefhet["VAF_2"] > 0.3
        ),
        table_vcf_mul_nonrefhet["GQ"] > 25
    )
]
# %%
############################
table_vcf["Length_Calc"] = table_vcf.apply(
    lambda row: len(row["ALT"]) - len(row["REF"]) if set(row["ALT"]) <= {'A','T','G','C'} else None, axis = 1
)

def get_snv_type_from_len(length):
    if length == 0:
        return "SNV"
    elif length < 0:
        return "DEL"
    else:
        return "INS"
def plot_counter(count_data, keys_sorted = None):
    if keys_sorted == None:
        keys_sorted = sorted(count_data.keys(), key = count_data.__getitem__, reverse = True)

    plt.bar(
        list(map(str, keys_sorted)),
        list(map(count_data.__getitem__, keys_sorted)),
        color = "gray"
    )
def convert_number_to_kmg(val):
    if val >= 1e9:
        val_annot = f"{val/1e9:.1f}G"
    elif val >= 1e6:
        val_annot = f"{val/1e6:.1f}M"
    elif val >= 1e3:
        val_annot = f"{val/1e3:.1f}K"
    else:
        val_annot = f"{val}"
    return val_annot

table_vcf["SNV_Type"] = table_vcf["Length_Calc"].apply(get_snv_type_from_len)

import numpy as np
from collections import Counter
from matplotlib import pyplot as plt
plt.rcParams["font.size"] = 15

count_bi_type = Counter(table_vcf["SNV_Type"])
plt.figure(figsize = (3.5, 4))
plot_counter(count_bi_type, keys_sorted = ["SNV", "INS", "DEL"])
for ind, key in enumerate(["SNV", "INS", "DEL"]):
    plt.annotate(convert_number_to_kmg(count_bi_type[key]), (ind, count_bi_type[key]), va = "bottom", ha = "center")
plt.ylim(top = max(count_bi_type.values())*1.1)
plt.xlim(-0.75, 2.75)
plt.show()
# %%
from collections import Counter
from matplotlib import pyplot as plt
import numpy as np
list_chrnames = list(map(lambda val: f"chr{val}", range(1, 23))) + ["chrX", "chrY"]

count_vcfmerfin_chrname = Counter(table_vcf["#CHROM"])
count_vcfmerfin_snv_chrname = Counter(table_vcf[table_vcf["SNV_Type"] == "SNV"]["#CHROM"])
count_vcfmerfin_ins_chrname = Counter(table_vcf[table_vcf["SNV_Type"] == "INS"]["#CHROM"])
count_vcfmerfin_del_chrname = Counter(table_vcf[table_vcf["SNV_Type"] == "DEL"]["#CHROM"])

plt.figure(figsize = (14, 3))
plt.rcParams['font.size'] = 20

arr_counts = np.zeros(len(list_chrnames))

plt.bar(list_chrnames, list(map(lambda key: count_vcfmerfin_ins_chrname[key], list_chrnames)), color = "cornflowerblue", bottom = arr_counts, label = "Insertion")

arr_counts += np.array(list(map(lambda key: count_vcfmerfin_ins_chrname[key], list_chrnames)))

plt.bar(list_chrnames, list(map(lambda key: count_vcfmerfin_del_chrname[key], list_chrnames)), color = "salmon", bottom = arr_counts, label = "Deletion")

arr_counts += np.array(list(map(lambda key: count_vcfmerfin_del_chrname[key], list_chrnames)))

plt.bar(list_chrnames, list(map(lambda key: count_vcfmerfin_snv_chrname[key], list_chrnames)), color = "limegreen", bottom = arr_counts, label = "SNV")


plt.xlabel("")
plt.ylabel("# SVs")
plt.xticks(rotation = 90)
plt.xlim(-0.5, len(list_chrnames)-0.5)
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1], loc = "center left", bbox_to_anchor = (1.01, 0.5))

for ind, chrname in enumerate(list_chrnames):
    plt.annotate(convert_number_to_kmg(count_vcfmerfin_chrname[chrname]), (ind, count_vcfmerfin_chrname[chrname]), va = "bottom", ha = "center")
plt.ylim(top = max(count_vcfmerfin_chrname.values())*1.2)
plt.show()
# %%
