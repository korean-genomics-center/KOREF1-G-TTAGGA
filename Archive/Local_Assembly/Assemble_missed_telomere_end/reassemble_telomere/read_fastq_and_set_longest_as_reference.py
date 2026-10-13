#%%
import os
from copy import deepcopy


def run(path_fastq, path_save_ref, path_save_others):
    dict_fastq = read_fastq(path_fastq)
    
    longest_read = sorted(
        dict_fastq.keys(), key = lambda name: len(dict_fastq[name]["seq"])
    )[-1]
    
    dict_fastq_wo_longest = deepcopy(dict_fastq)
    dict_seq_longest = dict_fastq_wo_longest.pop(longest_read)
    
    write_single_read_as_fasta(longest_read, dict_seq_longest, path_save_ref)
    write_fastq(path_save_others, dict_fastq_wo_longest)    

def read_fastq(path_fastq):
    dict_fastq = dict()
    
    with open(path_fastq, 'r') as fr:
        while 1:
            header = fr.readline()
            if not header:
                break
            
            header = header.strip('\n')[1:]
            sequence = fr.readline().strip('\n')
            _ = fr.readline()
            quality = fr.readline().strip('\n')
            
            dict_fastq[header] = {
                "seq":sequence,
                "qual":quality
            }
    return dict_fastq

def write_single_read_as_fasta(header, dict_seq, path_save):
    with open(path_save, 'w') as fw:
        fw.write(f">{header}\n")
        fw.write(f"{dict_seq['seq']}\n")

def write_fastq(path_write, dict_fastq):
    with open(path_write, 'w') as fw:
        for readid in dict_fastq.keys():
            fw.write(f"@{readid}\n")
            fw.write(f"{dict_fastq[readid]['seq']}\n")
            fw.write("+\n")
            fw.write(f"{dict_fastq[readid]['qual']}\n")
            
#%%
if __name__ == "__main__":
    for ind in range(7):
        path_fastq = f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/Cluster_Unaligned_ONT_Reads_by_Alignment/ONT_Reads.Telomeric_Read_Self_bothhap_a2h80idy95_ratio_60.mapontPcw5U11000000e0.idy95align1000hang100000a2h90.lenover_100000.alignment_ratio_90.readset_{ind}.fastq"
        
        path_save_ref = f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/Cluster_Unaligned_ONT_Reads_by_Alignment/ONT_Reads.Telomeric_Read_Self_bothhap_a2h80idy95_ratio_60.mapontPcw5U11000000e0.idy95align1000hang100000a2h90.lenover_100000.alignment_ratio_90.readset_{ind}.longest.fa"
        path_save_others = f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/Cluster_Unaligned_ONT_Reads_by_Alignment/ONT_Reads.Telomeric_Read_Self_bothhap_a2h80idy95_ratio_60.mapontPcw5U11000000e0.idy95align1000hang100000a2h90.lenover_100000.alignment_ratio_90.readset_{ind}.others.fastq"
    
        run(path_fastq, path_save_ref, path_save_others)