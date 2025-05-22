#%%
import pandas as pd
import os
from fasta_parser import read_fasta_as_dict, save_specific_seqnames_from_fasta_dict
#%%
hap = 2

# path_telo = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/Telomere_Detection/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.asm.dip.hap{hap}.p_ctg.tidk_search_CCCTAA.w_1000_telomeric_repeat_windows.tsv"
path_telo = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/Telomere_Detection/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg.tidk_search_CCCTAA.w_1000_telomeric_repeat_windows.tsv"
# path_telo = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Results/Telomere_Detection/CHM13.tidk_search_CCCTAA.w_1000_telomeric_repeat_windows.tsv"

# path_fasta = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.asm.dip.hap{hap}.p_ctg.fasta"
path_fasta = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg.fasta"

# path_save_fasta_prefix = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/Telomere_Detection/Fasta_Split_by_Telomere/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.asm.dip.hap{hap}.p_ctg"
path_save_fasta_prefix = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/Telomere_Detection/Fasta_Split_by_Telomere/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg"


os.makedirs(os.path.dirname(path_save_fasta_prefix), exist_ok = True)
# %%
table_telo = pd.read_csv(path_telo, sep = '\t')
# %%
def detect_telomere_end(table_telo, n_check_window = 5, telo_classify_limit = 100):
    dict_seqname_to_telo = dict()
    for seqname in table_telo["id"].unique():
        table_telo_seq = table_telo[table_telo["id"] == seqname].sort_values(by = ["window"]).reset_index(drop = True)
        
        table_telo_seq_front = table_telo_seq.loc[:n_check_window-1, :]
        table_telo_seq_end = table_telo_seq.loc[table_telo_seq.shape[0]-n_check_window:, :]

        sum_forward_repeat_front = table_telo_seq_front["forward_repeat_number"].sum()
        sum_reverse_repeat_front = table_telo_seq_front["reverse_repeat_number"].sum()

        sum_forward_repeat_end = table_telo_seq_end["forward_repeat_number"].sum()
        sum_reverse_repeat_end = table_telo_seq_end["reverse_repeat_number"].sum()

        dict_seqname_to_telo[seqname] = {
            "front" : {
                "forward" : sum_forward_repeat_front > telo_classify_limit,
                "reverse" : sum_reverse_repeat_front > telo_classify_limit
            },
            "end" : {
                "forward" : sum_forward_repeat_end > telo_classify_limit,
                "reverse" : sum_reverse_repeat_end > telo_classify_limit
            }
        }
    return dict_seqname_to_telo
#%%
dict_seqname_to_teloend = detect_telomere_end(table_telo)
#%%
list_seqnames_with_telomere = list(filter(lambda seqname: sum(dict_seqname_to_teloend[seqname]["front"].values())+sum(dict_seqname_to_teloend[seqname]["end"].values()) > 0, dict_seqname_to_teloend.keys()))
#%%
print(len(list_seqnames_with_telomere))
print(len(list(filter(lambda name: name not in list_seqnames_with_telomere, dict_seqname_to_teloend.keys()))))
#%%
from copy import deepcopy

# %%
list_seqnames_both_telo = list()
list_seqnames_single_telo = list()

for seqname in list_seqnames_with_telomere:
    front_fw = dict_seqname_to_teloend[seqname]["front"]["forward"]
    front_rv = dict_seqname_to_teloend[seqname]["front"]["reverse"]
    end_fw = dict_seqname_to_teloend[seqname]["end"]["forward"]
    end_rv = dict_seqname_to_teloend[seqname]["end"]["reverse"]

    assert (front_fw and end_fw) == False, seqname
    assert (front_rv and end_rv) == False, seqname
    assert (front_fw and front_rv) == False, seqname
    assert (end_fw and end_rv) == False, seqname

    is_front_telo = front_fw or front_rv
    is_end_telo = end_fw or end_rv

    if is_front_telo and is_end_telo:
        list_seqnames_both_telo.append(seqname)
    elif is_front_telo or is_end_telo:
        list_seqnames_single_telo.append(seqname)
    else:
        raise Exception(f"{seqname} doesnot contains telomere")
# %%
print(len(list_seqnames_both_telo))
print(len(list_seqnames_single_telo))

list_seqnames_no_telo = list(filter(lambda name: name not in list_seqnames_with_telomere, dict_seqname_to_teloend.keys()))
print(len(list_seqnames_no_telo))

# %%
#%%
from matplotlib import pyplot as plt
import math

list_seqnames_draw = list_seqnames_single_telo

n_ax_x = 3
n_ax_y = math.ceil(len(list_seqnames_draw)/n_ax_x)
size_ax_x = 6
size_ax_y = 3
linewidth = 3

ylim = (-190, 190)

plt.rcParams["font.size"] = 14
fig, ax = plt.subplots(n_ax_y, n_ax_x, figsize = (size_ax_x * n_ax_x, size_ax_y * n_ax_y))

axes = ax.flatten()

for ind, seqname in enumerate(list_seqnames_draw):
    table_telo_chr = table_telo[table_telo["id"] == seqname]

    axes[ind].plot(table_telo_chr["window"], table_telo_chr["forward_repeat_number"], color = "cornflowerblue", label = "forward (CCCTAA)", linewidth = linewidth)
    axes[ind].plot(table_telo_chr["window"], -1 * table_telo_chr["reverse_repeat_number"], color = "salmon", label = "reverse (TTAGGG)", linewidth = linewidth)

    if ind >= (n_ax_x *  (n_ax_y-1)):
        axes[ind].set_xlabel("Coordinates (bp)")
    if ind % n_ax_x == 0:
        axes[ind].set_ylabel("# of Repeats")
    axes[ind].set_title(seqname)
    
    if ind >= (n_ax_x * (n_ax_y-1)) and ind % n_ax_x == 0:
        axes[ind].legend(loc = "upper center", bbox_to_anchor = (0.5, -0.5))
    axes[ind].set_ylim(ylim)

plt.tight_layout()
plt.show()
# %%
dict_fasta = read_fasta_as_dict(path_fasta)

save_specific_seqnames_from_fasta_dict(
    dict_fasta,
    list_seqnames_both_telo,
    f"{path_save_fasta_prefix}.telo_both.fasta"
)
save_specific_seqnames_from_fasta_dict(
    dict_fasta,
    list_seqnames_single_telo,
    f"{path_save_fasta_prefix}.telo_single.fasta"
)
save_specific_seqnames_from_fasta_dict(
    dict_fasta,
    list_seqnames_no_telo,
    f"{path_save_fasta_prefix}.no_telo.fasta"
)
# %%
