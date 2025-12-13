#%%
import math

import pandas as pd
from joblib import Parallel, delayed


#%%
path_paf = "/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/Consensus_read/Minimap2_Align/Align_Minimap2.REF_Hap1.QRY_Racon_contig_unclustered.mapontPcw5r2ke0nolongjoin.idy95align1000hang1000000a2h0.paf"

path_save_aligned_ratio = "/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/Consensus_read/Minimap2_Align/Align_Minimap2.REF_Hap1.QRY_Racon_contig_unclustered.mapontPcw5r2ke0nolongjoin.idy95align1000hang1000000a2h0.alignment_ratio.txt"

n_jobs = 64

#%%
def read_paf(path_paf):
    table_paf = pd.read_csv(path_paf, sep = '\t', names = list(range(1, 13)), usecols = list(range(1, 13)), low_memory = False)
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

def collect_aligned_ratio_of_contigs_on_paf(table_paf_check, n_jobs):
    list_contignames_total = list(table_paf_check["qryname"].unique())
    list_of_list_contignames_for_parallel = split_values_by_n_batch(list_contignames_total, n_jobs)
    
    list_tables_for_parallel = list(map(lambda names: table_paf_check.loc[names, :], list_of_list_contignames_for_parallel))
    
    with Parallel(n_jobs = n_jobs) as parallel:
        list_dict_counts = parallel(
            delayed(collect_aligned_ratio_of_contigs_on_paf_batches)(
                table_check_batch
            ) for table_check_batch in list_tables_for_parallel
        )
    
    dict_contig_to_aligned_ratio = dict()
    for dict_batch_result in list_dict_counts:
        dict_contig_to_aligned_ratio.update(dict_batch_result)
    
    return dict_contig_to_aligned_ratio

def collect_aligned_ratio_of_contigs_on_paf_batches(table_paf_check):
    dict_contig_to_aligned_ratio = dict()

    for contigname in table_paf_check["qryname"].unique():
        table_paf_contig = table_paf_check.loc[[contigname], :].copy()
        
        contig_length = table_paf_contig["qrylen"].to_list()[0]
        
        dict_contig_to_aligned_ratio[contigname] = {
            "length": contig_length,
            "aligned_to_ref":dict(),
            "alignedratio_to_ref":dict()
        }
        
        for aligned_refname in table_paf_contig["refname"].unique():
            table_paf_contig_to_ref = table_paf_contig[table_paf_contig["refname"] == aligned_refname]
            
            list_aligned_positions = list()
            for _, row in table_paf_contig_to_ref.iterrows():
                positions_aligned = list(range(row["qry_alignstart"], row["qry_alignend"]))
                list_aligned_positions.extend(positions_aligned)
            
            set_aligned_positions = set(list_aligned_positions)
            
            length_aligned_regions = len(set_aligned_positions)
            
            dict_contig_to_aligned_ratio[contigname]["aligned_to_ref"][aligned_refname] = length_aligned_regions
            dict_contig_to_aligned_ratio[contigname]["alignedratio_to_ref"][aligned_refname] = length_aligned_regions / contig_length
    return dict_contig_to_aligned_ratio

def split_values_by_n_batch(list_values, n_split):
    len_values = len(list_values)

    n_values_per_split = math.floor(len_values / n_split)
    n_leftovers = len_values - n_values_per_split * n_split
    
    list_values_split = list()
    
    start = 0
    for ind in range(n_split):
        n_values_this_ind = n_values_per_split + int(ind < n_leftovers)
        list_values_split.append(list_values[start:start+n_values_this_ind])
        
        start += n_values_this_ind
    return list_values_split
    
    
#%%
table_paf = read_paf(path_paf)
#%%
table_paf_indexed_by_contig = table_paf.set_index("qryname", drop = False)
# %%
dict_contig_to_aligned_ratio = collect_aligned_ratio_of_contigs_on_paf(table_paf_indexed_by_contig, n_jobs)
# %%
with open(path_save_aligned_ratio, 'w') as fw:
    for qryname in dict_contig_to_aligned_ratio.keys():
        lenqry = dict_contig_to_aligned_ratio[qryname]["length"]
        
        for refname in dict_contig_to_aligned_ratio[qryname]["aligned_to_ref"].keys():
            refalign_len = dict_contig_to_aligned_ratio[qryname]["aligned_to_ref"][refname]
            refalign_ratio = dict_contig_to_aligned_ratio[qryname]["alignedratio_to_ref"][refname]
            
            fw.write(f"{qryname}\t{lenqry}\t{refname}\t{refalign_len}\t{refalign_ratio}\n")
# %%
