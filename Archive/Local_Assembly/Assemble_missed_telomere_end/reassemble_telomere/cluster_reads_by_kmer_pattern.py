#%%
import pandas as pd
import numpy as np
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt
import seaborn as sns
from sklearn.cluster import OPTICS
from sklearn.cluster import DBSCAN
from collections import Counter
#%%
path_count = "/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/Classify_ONT_reads/KoGoral.Merged.ONT_reads.Adapt_Filt.Edge_Trim_Front_500_Tail_500.Filter_L50kb_QV15.telomere_count.highrepeat_500.readid.unaligned_reads.ratio_65.k17.min_12.max_42.count"

path_save_tsne = "/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/Classify_ONT_reads/KoGoral.Merged.ONT_reads.Adapt_Filt.Edge_Trim_Front_500_Tail_500.Filter_L50kb_QV15.telomere_count.highrepeat_500.readid.unaligned_reads.ratio_65.k17.min_12.max_42.count.tsne"

ncomp = 2
perplexity = 5
#%%
table_count = pd.read_csv(path_count, sep = '\t', names = ["readid", "readlen", "kmer", "kmer_count"])

# %%
table_count_matrix = table_count.pivot_table(index = "readid", columns = "kmer", values = "kmer_count")
#%%
table_existence_matrix = table_count_matrix.applymap(lambda val: int(val > 0))
# %%
table_tsne = TSNE(n_components=ncomp, learning_rate="auto", init="pca", perplexity=perplexity).fit_transform(table_existence_matrix)

#%%
table_tsne_df = pd.DataFrame(table_tsne).rename(columns = {0:"tSNE1", 1:"tSNE2"})
table_tsne_df["readid"] = table_existence_matrix.index
# %%

optics_clustering = OPTICS(min_samples=5).fit(table_tsne_df[["tSNE1", "tSNE2"]])
# hdbscan_clustering = DBSCAN(min_samples = 5, eps = 5).fit(table_tsne_df[["tSNE1", "tSNE2"]])
# hdbscan_clustering.labels_
# %%
table_tsne_df["Cluster"] = optics_clustering.labels_
#%%

# fig, ax = plt.subplots(2, 2, figsize = (7, 7))

# palette = sns.color_palette("tab10")
# dict_palette = dict(zip(sorted(table_tsne_df["Cluster"].unique()), palette))


# for ts1 in range(1, 3):
#     for ts2 in range(2, 4):
        
#         if ts1 == ts2: continue
        
#         sns.scatterplot(
#             data = table_tsne_df,
#             x = f"tSNE{ts1}",
#             y = f"tSNE{ts2}",
#             hue = "Cluster",
#             palette=dict_palette,
#             ax = ax[ts2-2, ts1-1]
#         )
# %%
sns.scatterplot(
    data = table_tsne_df,
    x = f"tSNE1",
    y = f"tSNE2",
    hue = "Cluster",
    palette="tab10"
)
# %%
table_tsne_df.to_csv(path_save_tsne, sep = '\t', index = False)
# %%
