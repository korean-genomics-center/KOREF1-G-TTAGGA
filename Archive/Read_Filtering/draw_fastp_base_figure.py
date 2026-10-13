#%%
import json
from matplotlib import pyplot as plt
import seaborn as sns
#%%
sampleid = "PATERNAL"
runs = list(range(1, 9))

dict_quality = {
    "read1":list(),
    "read2":list()
}
dict_gc = {
    "read1":list(),
    "read2":list()
}

total_base_bef = 0
total_base_aft = 0

list_read1_length = list()
list_read2_length = list()

for run in runs:
    path_report = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/{sampleid}/Trim_FQ/{sampleid}_L{run}_Trimmed_F15_T5.noNbase.report.json"
    with open(path_report, 'r') as fr:
        dict_report = json.load(fr)
        
        list_read1_quality = dict_report["read1_after_filtering"]["quality_curves"]["mean"]
        list_read2_quality = dict_report["read2_after_filtering"]["quality_curves"]["mean"]
        list_read1_gc = dict_report["read1_after_filtering"]["content_curves"]["GC"]
        list_read2_gc = dict_report["read2_after_filtering"]["content_curves"]["GC"]
        
        total_base_bef += dict_report["summary"]["before_filtering"]["total_bases"]
        total_base_aft += dict_report["summary"]["after_filtering"]["total_bases"]
        
        list_read1_length.append(dict_report["summary"]["after_filtering"]["read1_mean_length"])
        list_read2_length.append(dict_report["summary"]["after_filtering"]["read2_mean_length"])
        
        dict_quality["read1"].append(list_read1_quality)
        dict_quality["read2"].append(list_read2_quality)
        dict_gc["read1"].append(list_read1_gc)
        dict_gc["read2"].append(list_read2_gc)
# %%
plt.figure(figsize = (10, 4))

dict_show = dict_quality
read_show = "read2"
for ind, data in enumerate(dict_show[read_show]):
    plt.plot(range(1, len(data)+1), data, color = sns.color_palette("hls", len(runs))[ind])

plt.xlabel("Coordinate (bp)")
plt.ylabel("Quality")
plt.show()
# %%
plt.figure(figsize = (10, 4))

dict_show = dict_gc
read_show = "read2"
for ind, data in enumerate(dict_show[read_show]):
    plt.plot(range(1, len(data)+1), data, color = sns.color_palette("hls", len(runs))[ind])
    
    if sum(list(map(lambda val: val < 0.36, data))) > 0:
        print(ind)

plt.xlabel("Coordinate (bp)")
plt.ylabel("GC Ratio")
plt.ylim(0, 0.5)
plt.show()
# %%
