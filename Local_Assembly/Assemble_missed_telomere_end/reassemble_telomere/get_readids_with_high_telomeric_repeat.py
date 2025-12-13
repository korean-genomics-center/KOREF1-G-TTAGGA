#%%
import pandas as pd
import numpy as np

path_repeat_count = "/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/Collect_ONT_with_telomere/KoGoral.Merged.ONT_reads.Adapt_Filt.Edge_Trim_Front_500_Tail_500.Filter_L50kb_QV15.telomere_count.txt"

highrepeat_thresh = 500
path_save_filtered_readids = f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/Collect_ONT_with_telomere/KoGoral.Merged.ONT_reads.Adapt_Filt.Edge_Trim_Front_500_Tail_500.Filter_L50kb_QV15.telomere_count.highrepeat_{highrepeat_thresh}.readid.txt"

#%%
table_repeat_count = pd.read_csv(path_repeat_count, sep = '\t', names = ["readid", "length", "fw_count", "rv_count"])
# %%
table_repeat_count["fw_highrepeat"] = table_repeat_count["fw_count"] > highrepeat_thresh
table_repeat_count["rv_highrepeat"] = table_repeat_count["rv_count"] > highrepeat_thresh
# %%
table_repeat_count_highrepeat = table_repeat_count[np.logical_or(table_repeat_count["fw_highrepeat"], table_repeat_count["rv_highrepeat"])]
# %%
with open(path_save_filtered_readids, 'w') as fw:
    for readid in table_repeat_count_highrepeat["readid"]:
        fw.write(f"{readid}\n")
# %%
