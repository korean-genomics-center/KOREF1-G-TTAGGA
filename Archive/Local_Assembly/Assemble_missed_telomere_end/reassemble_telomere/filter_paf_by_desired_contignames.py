#%%
import pandas as pd


path_teloinfo_format = "/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Initial_Assembly/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_Detection/Fasta_Split_by_Telomere/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.asm.hic.hap{hap}.p_ctg.telo_single.telo_info.tsv"


path_paf_format = "/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/Consensus_read/Minimap2_Align/Align_Minimap2.REF_Hap{hap}.QRY_Racon_contig_ind{ind_contig}.mapontPcw5r2ke0nolongjoin.idy95align1000hang1000000a2h0.paf"


for hap in [1,2]:
    for ind_contig in [0,1,2,3,4,6]:
        path_paf_check = path_paf_format.format(hap = hap, ind_contig = ind_contig)
        
        set_contigs_with_telo = set(pd.read_csv(path_teloinfo_format.format(hap = hap), sep = '\t')["Scaffold"].to_list())
        
        path_paf_save_filt = path_paf_check.replace(".paf", ".telocontig.paf")
        
        with open(path_paf_check, 'r') as fr, open(path_paf_save_filt, 'w') as fw:
            for line in fr:
                contig_aligned = line.split('\t')[5]
                if contig_aligned in set_contigs_with_telo:
                    fw.write(line)