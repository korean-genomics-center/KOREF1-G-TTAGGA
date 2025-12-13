#%%
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
#%%

list_path_alignment_ratio = [
    "/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/Collect_ONT_with_telomere/Minimap2_Align/Align_Minimap2.REF_Hap1.QRY_TelomericONT.mapontPcw5r2ke0nolongjoin.idy95align1000hang100000a2h80.alignment_ratio.txt",
    "/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/Collect_ONT_with_telomere/Minimap2_Align/Align_Minimap2.REF_Hap2.QRY_TelomericONT.mapontPcw5r2ke0nolongjoin.idy95align1000hang100000a2h80.alignment_ratio.txt"
]

alignratio_high = 0.60

path_read = "/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/Collect_ONT_with_telomere/KoGoral.Merged.ONT_reads.Adapt_Filt.Edge_Trim_Front_500_Tail_500.Filter_L50kb_QV15.telomere_count.highrepeat_500.readid.txt"
path_save_readid = f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/Collect_ONT_with_telomere/KoGoral.Merged.ONT_reads.Adapt_Filt.Edge_Trim_Front_500_Tail_500.Filter_L50kb_QV15.telomere_count.highrepeat_500.readid.unaligned_reads.bothhap.a2h80idy95.ratio_{int(alignratio_high*100)}.txt"

#%%

table_alignment_ratio = pd.concat(
    list(map(lambda path: pd.read_csv(path, sep = '\t', names = ["readid", "readlen", "contigid", "alignlen", "alignratio"]), list_path_alignment_ratio))
)
# %%
table_alignment_ratio_onlyhighest = table_alignment_ratio.sort_values(by = ["alignratio"], ascending=False).drop_duplicates(subset = ["readid"])
# %%
plt.rcParams["font.size"] = 14
plt.hist(table_alignment_ratio_onlyhighest["alignratio"], bins = np.linspace(0, 1, 101), color = "gray")
plt.xlabel("Aligned Ratio")
plt.ylabel("# Reads")
# %%
readids_high_alignratio = table_alignment_ratio_onlyhighest[table_alignment_ratio_onlyhighest["alignratio"] > alignratio_high]["readid"].to_list()
# %%
with open(path_read, 'r') as fr, open(path_save_readid, 'w') as fw:
    for line in fr:
        readid = line.strip('\n')
        if readid in readids_high_alignratio:
            continue
        else:
            fw.write(readid + '\n')
# %%
