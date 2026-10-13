#%%
import math

from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

# %%
path_flagger_result = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Flagger_Evaluation/Flagger/ONT_10kb_QV15/prediction_summary_final.tsv"

table_flagger_result = pd.read_csv(path_flagger_result, sep = '\t')
# %%
table_flagger_result_final = table_flagger_result[
    np.logical_and(table_flagger_result["Metric_Type"] == "base_level", table_flagger_result["Category_Type"] == "annotation")
]
# %%
table_flagger_result_final_perc = table_flagger_result_final[table_flagger_result_final["Entry_Type"] == "percentage"]
# %%
cols_show = [
    "WHOLE_GENOME",
    "CNTR_or_SD",
    "Non_CNTR_and_Non_SD",
    "AUTOSOME",
    "AUTOSOME_Non_CNTR_and_Non_SD"
]

def get_renamed_val(val):
    list_val = val.split('_')

    list_rename_val = list()
    for val in list_val:
        if val.lower() in ["sd"]:
            val = val
        elif val.lower() in ["or", "and"]:
            val = val.lower()
        else:
            val = val[0].upper() + val[1:].lower()
        list_rename_val.append(
            val
        )
    return '_'.join(list_rename_val)
#%%
dict_name_to_hap = dict(zip(table_flagger_result_final_perc["Category_Name"], table_flagger_result_final_perc["Hap"]))
dict_name_to_err = dict(zip(table_flagger_result_final_perc["Category_Name"], table_flagger_result_final_perc["Err"]))
dict_name_to_dup = dict(zip(table_flagger_result_final_perc["Category_Name"], table_flagger_result_final_perc["Dup"]))
dict_name_to_col = dict(zip(table_flagger_result_final_perc["Category_Name"], table_flagger_result_final_perc["Col"]))
# %%
plt.rcParams["font.size"] = 14
plt.rcParams["font.family"] = "sans-serif"
# %%
plt.figure(figsize = (4, 3))
plt.bar(
    list(map(get_renamed_val, cols_show)),
    list(map(dict_name_to_hap.__getitem__, cols_show)),
    color = "gray"
)
for ind, name in enumerate(cols_show):
    plt.annotate(dict_name_to_hap[name], (ind, dict_name_to_hap[name]), va = "bottom", ha = "center")
plt.xticks(rotation = -45, rotation_mode="anchor", ha = "left", va = "top")
plt.ylim(top = 110)
plt.ylabel("Correct Regions (%)")
# %%
# %%
plt.figure(figsize = (4, 3))

dict_issue_to_color = {
    "Dup":(250/255,104/255,0/255),
    "Col":(170/255,0/255,255/255),
    "Err":(162/255,0/255,37/255)
}

bottom = np.zeros(len(cols_show))
for issue in ["Err", "Col", "Dup"]:
    dict_name_to_issue = dict(zip(table_flagger_result_final_perc["Category_Name"], table_flagger_result_final_perc[issue]))

    values_issue = np.array(list(map(dict_name_to_issue.__getitem__, cols_show)))
    plt.bar(
        list(map(get_renamed_val, cols_show)),
        values_issue,
        color = dict_issue_to_color[issue],
        bottom = bottom,
        label = issue
    )
    bottom += values_issue
plt.xticks(rotation = -45, rotation_mode="anchor", ha = "left", va = "top")
plt.ylim(top = 4)
plt.ylabel("Issued Regions (%)")
plt.legend(loc = "center left", bbox_to_anchor = (1.01, 0.5))
# %%
