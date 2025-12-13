#%%
import os, sys


sys.path.append("/BiO/Research/Project2/Endangered_Goral/Resources/Scripts/KOGIC_KOREF_AGGA/Gap_Filling/Manual_filling")
from check_contig_overlap_elongation import *


#%%
path_paf = "//BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/From_ONT/RawReads_Aligned_ONT/from_round1/Minimap2_Align/Align_Minimap2.Ref_h1tg000026l.R1.QRY_Raw_ONT.mapontPcw10U101000000.endaligned1000.tmp.idy95align1000hang1000a2h95.paf"


dict_contigs_collect = {
    # "h1tg000011l_elongate_1":"end",
    # "h1tg000012l_elongate_1":"end",
    "h1tg000026l_elongate_1":"end",
}
# %%
table_paf = read_paf(path_paf)
#%%
table_paf_filt = filter_paf_of_end_aligned_query(table_paf, min_idy = 0.99, min_align = 10000, max_hang = 100, min_align2hang = 0.99)
# %%
dict_contig_to_elong_alignment = dict()
dict_best_read = dict()

for contig, direction in dict_contigs_collect.items():
    table_paf_filt_contig = table_paf_filt[table_paf_filt["refname"] == contig]
    
    if direction == "front":
        in_range = table_paf_filt_contig["ref_alignstart"] <= 100
    else:
        reflen = table_paf_filt_contig["reflen"].to_list()[0]
        in_range = table_paf_filt_contig["ref_alignend"] >= (reflen-100)
    
    table_paf_filt_contig_for_elong = table_paf_filt_contig[in_range]
    
    dict_contig_to_elong_alignment[contig] = table_paf_filt_contig_for_elong.copy()
    print(len(table_paf_filt_contig_for_elong["qryname"].unique()))
    
    table_paf_filt_contig_for_elong["elongation"] = table_paf_filt_contig_for_elong["qrylen"] - table_paf_filt_contig_for_elong["align"]
    
    table_paf_filt_contig_for_elong_sort = table_paf_filt_contig_for_elong.sort_values(by = ["elongation"], ascending = False)
    
    best_alignment = table_paf_filt_contig_for_elong_sort.iloc[0, :]
    dict_best_read[contig] = best_alignment

#%%

    
# %%
for contig in dict_contigs_collect.keys():
    with open(f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/From_ONT/RawReads_Aligned_ONT/from_contig/{contig}_end.ont_readid.txt", 'w') as fw:
        for readid in dict_contig_to_elong_alignment[contig]["qryname"].unique():
            fw.write(f"{readid}\n")
# %%
##################################
table_paf_filt["elongation"] = table_paf_filt["qrylen"] - table_paf_filt["align"]
table_paf_filt[table_paf_filt["refname"] == "h1tg000026l"]["elongation"].hist()
# %%

# %%
for contig in dict_contigs_collect.keys():
    with open(f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/From_ONT/RawReads_Aligned_ONT/from_contig/{contig}_end.single_best.ont_readid.txt", 'w') as fw:
        fw.write(f"{dict_best_read[contig]['qryname']}\n")
# %%
