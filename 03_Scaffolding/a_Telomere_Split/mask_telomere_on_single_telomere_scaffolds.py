#%%
import pandas as pd
import os, sys
from pathlib import Path

sys.path.append(str(Path(__file__).parents[2] / "common"))

from fasta_parser import read_fasta_as_dict, save_dict_fasta_into_file
# %%
window_check_start = 5
telo_classify_limit = 100
telo_elongation_limit = 20

hap = 2

# path_telo = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/Telomere_Detection/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.asm.dip.hap{hap}.p_ctg.tidk_search_CCCTAA.w_1000_telomeric_repeat_windows.tsv"
# path_fasta_single_telo = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/Telomere_Detection/Fasta_Split_by_Telomere/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.asm.dip.hap{hap}.p_ctg.telo_single.fasta"
# path_save_prefix = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/Telomere_Detection/Fasta_Split_by_Telomere/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.asm.dip.hap{hap}.p_ctg.telo_single"

path_telo = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/Telomere_Detection/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg.tidk_search_CCCTAA.w_1000_telomeric_repeat_windows.tsv"
path_fasta_single_telo = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/Telomere_Detection/Fasta_Split_by_Telomere/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg.telo_single.fasta"
path_save_prefix = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/Telomere_Detection/Fasta_Split_by_Telomere/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg.telo_single"

os.makedirs(os.path.dirname(path_save_prefix), exist_ok = True)
#%%
table_telo = pd.read_csv(path_telo, sep = '\t')
dict_fasta = read_fasta_as_dict(path_fasta_single_telo)

#%%
table_telo_for_check = table_telo[table_telo["id"].isin(list(dict_fasta.keys()))]
#%%
def detect_telomere_end(table_telo, n_check_window = 5, telo_classify_limit = 100):
    dict_seqname_to_telo = dict()
    for seqname in table_telo["id"].unique():
        table_telo_seq = table_telo[table_telo["id"] == seqname].sort_values(by = ["window"]).reset_index(drop = True)
        
        table_telo_seq_front = table_telo_seq.loc[:n_check_window-1, :]
        table_telo_seq_end = table_telo_seq.loc[table_telo_seq.shape[0]-n_check_window:, :]

        print(table_telo_seq_front.shape)
        print(table_telo_seq_end.shape)

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
dict_seqname_to_teloend = detect_telomere_end(table_telo_for_check, n_check_window = window_check_start, telo_classify_limit = telo_classify_limit)
# %%
def get_telomere_sequence_window_ranges(table_telo, dict_seqname_to_teloend, window_check_start = 5, telo_elongation_limit = 20, telo_classify_limit = 100):
    dict_seqname_to_telorange = dict()
    for seqname in dict_seqname_to_teloend.keys():
        table_telo_seq = table_telo[table_telo["id"] == seqname].sort_values(by = ["window"]).reset_index(drop = True)

        is_telo_front = sum(dict_seqname_to_teloend[seqname]["front"].values()) == 1
        is_telo_end = sum(dict_seqname_to_teloend[seqname]["end"].values()) == 1
        assert is_telo_front or is_telo_end, seqname

        is_telo_forward = dict_seqname_to_teloend[seqname]["front"]["forward"] if is_telo_front else dict_seqname_to_teloend[seqname]["end"]["forward"]
        telo_count_colname = "forward_repeat_number" if is_telo_forward else "reverse_repeat_number"

        if is_telo_end:
            table_telo_seq = table_telo_seq.sort_index(ascending = False)
        
        n_elongation = 0
        assert table_telo_seq.iloc[:window_check_start, :][telo_count_colname].sum() >= telo_classify_limit, seqname
        while 1:
            row = table_telo_seq.iloc[window_check_start+n_elongation, :]
            if row[telo_count_colname] >= telo_elongation_limit:
                n_elongation += 1
            else:
                n_elongation -= 1
                break
    

        if is_telo_front:
            window_elongation_stopped = table_telo_seq.iloc[window_check_start+n_elongation, :]["window"]
            dict_seqname_to_telorange[seqname] = (0, window_elongation_stopped)
        else:
            window_before_elongation_stopped = table_telo_seq.iloc[window_check_start+n_elongation+1, :]["window"]
            last_window = table_telo_seq.iloc[0, :]["window"]
            dict_seqname_to_telorange[seqname] = (window_before_elongation_stopped, last_window)
    return dict_seqname_to_telorange
# %%
dict_seqname_to_telorange = get_telomere_sequence_window_ranges(table_telo_for_check, dict_seqname_to_teloend, window_check_start = window_check_start, telo_elongation_limit = telo_elongation_limit, telo_classify_limit = telo_classify_limit)

path_save_masked_fasta = f"{path_save_prefix}.telo_masked.fasta"
path_save_telo_fasta = f"{path_save_prefix}.masked_telo_only.fasta"

path_save_teloend_table = f"{path_save_prefix}.telo_info.tsv"


dict_seqname_to_masked_fasta = dict()
dict_seqname_to_telo_fasta = dict()

for seqname in dict_seqname_to_telorange.keys():
    rangestart, rangeend = dict_seqname_to_telorange[seqname]
    telo_seq = dict_fasta[seqname][rangestart:rangeend]

    seq_masked = dict_fasta[seqname].replace(telo_seq, 'N' * len(telo_seq))

    dict_seqname_to_masked_fasta[f"{seqname}_telomasked"] = seq_masked
    dict_seqname_to_telo_fasta[f"{seqname}_teloonly_{rangestart}_{rangeend}"] = telo_seq
#%%
table_teloend = pd.DataFrame({
    "Scaffold" : list(dict_seqname_to_teloend.keys())
})
table_teloend["Telomere_End"] = table_teloend["Scaffold"].apply(lambda seqname: "front" if sum(dict_seqname_to_teloend[seqname]["front"].values()) == 1 else "end")
table_teloend["Telomere_Sequence_Direction"] = table_teloend.apply(lambda row: "forward" if dict_seqname_to_teloend[row["Scaffold"]][row["Telomere_End"]]["forward"] else "reverse", axis = 1)

table_teloend["Scaffold_Length"] = table_teloend["Scaffold"].apply(lambda seqname: len(dict_fasta[seqname]))
table_teloend["Telomere_Masking_Start"] = table_teloend["Scaffold"].apply(lambda seqname: dict_seqname_to_telorange[seqname][0])
table_teloend["Telomere_Masking_End"] = table_teloend["Scaffold"].apply(lambda seqname: dict_seqname_to_telorange[seqname][1])

# %%
table_teloend.to_csv(path_save_teloend_table, sep = '\t', index = False)
save_dict_fasta_into_file(dict_seqname_to_masked_fasta, path_save_masked_fasta)
save_dict_fasta_into_file(dict_seqname_to_telo_fasta, path_save_telo_fasta)
# %%