#%%
import pandas as pd
import numpy as np
import networkx as nx
from matplotlib import pyplot as plt
import seaborn as sns
from networkx.algorithms import community
#%%
path_aligned_ratio = "/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/Cluster_Unaligned_ONT_Reads_by_Alignment/Minimap2_Align/Align_Minimap2.Telomeric_Read_Self_bothhap_a2h80idy95_ratio_60.mapontPcw5U11000000e0.idy95align1000hang100000a2h90.lenover_100000.alignment_ratio.txt"

path_save_readset_format = "/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/Cluster_Unaligned_ONT_Reads_by_Alignment/ONT_Reads.Telomeric_Read_Self_bothhap_a2h80idy95_ratio_60.mapontPcw5U11000000e0.idy95align1000hang100000a2h90.lenover_100000.alignment_ratio_90.readset_{ind}.txt"

path_save_leftovers = "/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/Cluster_Unaligned_ONT_Reads_by_Alignment/ONT_Reads.Telomeric_Read_Self_bothhap_a2h80idy95_ratio_60.mapontPcw5U11000000e0.idy95align1000hang100000a2h90.lenover_100000.alignment_ratio_90.noncluster.txt"

thresh_alignhigh = 0.9
#%%
table_aligned_ratio = pd.read_csv(path_aligned_ratio, sep = '\t', names = ["readid_q", "readlen_q", "readid_r", "alignlen", "alignratio"])

# %%
dict_read_aligned_pairwise = dict()

for _, row in table_aligned_ratio.iterrows():
    if row["readid_q"] == row["readid_r"]: continue
    readset = '__'.join(sorted([row["readid_q"], row["readid_r"]]))
    
    if dict_read_aligned_pairwise.get(readset) == None or dict_read_aligned_pairwise.get(readset) < row["alignratio"]:
        dict_read_aligned_pairwise[readset] = row["alignratio"]
# %%
dict_read_aligned_pairwise_filt = dict()

for readset, alignratio in dict_read_aligned_pairwise.items():
    if alignratio > thresh_alignhigh:
        dict_read_aligned_pairwise_filt[readset] = alignratio
# %%
vertex = set(table_aligned_ratio["readid_q"]) | set(table_aligned_ratio["readid_r"])

edge = list()
for readset in dict_read_aligned_pairwise_filt.keys():
    read1, read2 = readset.split("__")
    edge.append((read1, read2))


Graph = nx.Graph()
Graph.add_nodes_from(vertex)
Graph.add_edges_from(edge)
# %%
communities = list(community.greedy_modularity_communities(Graph))

node_to_community = {}
for i, comm in enumerate(communities):
    for node in comm:
        node_to_community[node] = i

pos = nx.spring_layout(Graph, seed=0)

colors = [node_to_community[node] for node in Graph.nodes()]

plt.figure(figsize=(5, 4))
nx.draw_networkx_nodes(Graph, pos, node_color=colors, cmap=plt.cm.jet, node_size=10)
nx.draw_networkx_edges(Graph, pos, alpha=0.5)

for comind in set(node_to_community.values()):
    list_reads_in_thiscom = list(filter(lambda name: node_to_community[name] == comind, node_to_community.keys()))
    
    list_xpos = list(map(lambda name: pos[name][0], list_reads_in_thiscom))
    list_ypos = list(map(lambda name: pos[name][1], list_reads_in_thiscom))
    
    plt.annotate(f"Ind{comind} (N={len(list_reads_in_thiscom)})", (np.mean(list_xpos), np.mean(list_ypos)), ha = "center", va = "bottom")

# nx.draw_networkx_labels(Graph, pos, font_size=0)
plt.show()
#%%
plt.figure(figsize=(10, 8))
nx.draw_networkx_nodes(Graph, pos, node_color=colors, cmap=plt.cm.jet, node_size=50)
nx.draw_networkx_edges(Graph, pos, alpha=0.5)
# nx.draw_networkx_labels(Graph, pos, font_size=0)
plt.show()

#%%

write_leftovers = open(path_save_leftovers, 'w')

ind_readset = 0
for comm in communities:
    if len(comm) >= 3:
        with open(path_save_readset_format.format(ind = ind_readset), 'w') as fw:
            for node in comm:
                fw.write(f"{node}\n")
        ind_readset += 1
    else:
        for node in comm:
            write_leftovers.write(f"{node}\n")
write_leftovers.close()
# %%
#################################
def read_paf(path_paf):
    table_paf = pd.read_csv(path_paf, sep = '\t', names = list(range(1, 13)), usecols = list(range(1, 13)))
    table_paf = table_paf.rename(
        columns = {
            1:"qryname",
            2:"qrylen",
            3:"qry_alignstart",
            4:"qry_alignend",
            5:"strand",
            6:"refname",
            7:"reflen",
            8:"ref_alignstart",
            9:"ref_alignend",
            10:"match",
            11:"align",
            12:"quality"
        }
    )
    return table_paf

table_paf = read_paf("/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/Cluster_Unaligned_ONT_Reads_by_Alignment/Minimap2_Align/Align_Minimap2.Telomeric_Read_Self_bothhap_a2h80idy95_ratio_60.mapontPcw5U11000000e0.idy95align1000hang100000a2h70.lenover_75000.alignover_75000.paf")
# %%
table_paf_wo_selfalign = table_paf[table_paf["qryname"] != table_paf["refname"]]

vertex = set(table_paf_wo_selfalign["qryname"]) | set(table_paf_wo_selfalign["refname"])

edge = list()
for _, row in table_paf_wo_selfalign.iterrows():
    if row["qryname"] == row["refname"]: continue
    vertex_set = sorted([row["qryname"], row["refname"]])
    if vertex_set in edge:
        continue
    else:
        edge.append(vertex_set)
# %%
Graph = nx.Graph()
Graph.add_nodes_from(vertex)
Graph.add_edges_from(edge)

communities = list(community.greedy_modularity_communities(Graph))

node_to_community = {}
for i, comm in enumerate(communities):
    for node in comm:
        node_to_community[node] = i

pos = nx.spring_layout(Graph, seed=0)

colors = [node_to_community[node] for node in Graph.nodes()]

plt.figure(figsize=(5, 4))
nx.draw_networkx_nodes(Graph, pos, node_color=colors, cmap=plt.cm.jet, node_size=10)
nx.draw_networkx_edges(Graph, pos, alpha=0.5)

for comind in set(node_to_community.values()):
    list_reads_in_thiscom = list(filter(lambda name: node_to_community[name] == comind, node_to_community.keys()))
    
    list_xpos = list(map(lambda name: pos[name][0], list_reads_in_thiscom))
    list_ypos = list(map(lambda name: pos[name][1], list_reads_in_thiscom))
    
    plt.annotate(f"Ind{comind} (N={len(list_reads_in_thiscom)})", (np.mean(list_xpos), np.mean(list_ypos)), ha = "center", va = "bottom")

# nx.draw_networkx_labels(Graph, pos, font_size=0)
plt.show()
# %%
from collections import Counter

Counter(node_to_community.values())

# %%
plt.figure(figsize=(10, 8))
nx.draw_networkx_nodes(Graph, pos, node_color=colors, cmap=plt.cm.jet, node_size=50)
nx.draw_networkx_edges(Graph, pos, alpha=0.5)
# nx.draw_networkx_labels(Graph, pos, font_size=0)
plt.show()
# %%
