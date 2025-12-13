#%%
import os, re, subprocess

import numpy as np
from joblib import Parallel, delayed


def run(path_file, path_save, n_jobs, known_repeat = "CCCTAA"):
    list_offsets_start_read = get_offsets_to_read_file_parallel(path_file, n_jobs)
    list_offsets_end_read = list_offsets_start_read[1:] + [list_offsets_start_read[-1]*10]
    
    list_path_save_tmp = list(map(lambda ind: f"{path_save}.tmp.{ind}", range(n_jobs)))
    with Parallel(n_jobs = n_jobs) as parallel:
        parallel(
            delayed(count_repeat_from_fasta_with_offsets)(
                path_file,
                list_path_save_tmp[ind_job],
                list_offsets_start_read[ind_job],
                list_offsets_end_read[ind_job],
                known_repeat
            ) for ind_job in range(n_jobs)
        )
    
    subprocess.run(f"cat {path_save}.tmp* > {path_save} && rm {path_save}.tmp*", shell = True)
    

def get_offsets_to_read_file_parallel(path_file, n_split):
    filesize = os.path.getsize(path_file)
    
    list_points_for_split = np.linspace(0, filesize, n_split+1)[:-1]
    
    list_offset_start_read = [0] + list(map(lambda point: adjust_startpoint_for_readstart(path_file, int(point)), list_points_for_split[1:]))
    return list_offset_start_read
    
def adjust_startpoint_for_readstart(path_file, offset):
    with open(path_file, 'r') as fr:
        fr.seek(offset)
        tmp = fr.readline()
        curr = fr.tell()
        line = fr.readline()
        next = fr.tell()
        
    if line.startswith('>'):
        return curr
    else:
        return next
    
def get_reverse_complement(sequence):
    dict_reverse = {
        'A':'T',
        'T':'A',
        'G':'C',
        'C':'G'
    }
    
    reverse_sequence = list(map(dict_reverse.__getitem__, sequence))
    rc_sequence = ''.join(reverse_sequence[::-1])
    return rc_sequence

def count_repeat_from_fasta_with_offsets(path_file, path_save_tmp, offset_start, offset_end, seq):
    rc_seq = get_reverse_complement(seq)
    
    with open(path_file, 'r') as fr, open(path_save_tmp, 'w') as fw:
        fr.seek(offset_start)
        
        while 1:
            header = fr.readline()
            if not header:
                break
            header = header.strip('\n')[1:].split()[0]
            sequence = fr.readline().strip('\n')
            
            count_fw_seq = sequence.count(seq)
            count_rv_seq = sequence.count(rc_seq)
            
            fw.write(f"{header}\t{len(sequence)}\t{count_fw_seq}\t{count_rv_seq}\n")
            
            if fr.tell() >= offset_end:
                break
# %%
if __name__ == "__main__":
    path_file = "/BiO/Research/Project2/Endangered_Goral/Results/Input_Data_Filtering/ONT/KoGoral/Merged_Fastq/KoGoral.Merged.ONT_reads.Adapt_Filt.Edge_Trim_Front_500_Tail_500.Filter_L50kb_QV15.fasta"
    path_save = "/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/Collect_ONT_with_telomere/KoGoral.Merged.ONT_reads.Adapt_Filt.Edge_Trim_Front_500_Tail_500.Filter_L50kb_QV15.telomere_count.txt"
    
    n_jobs = 64
    
    run(path_file, path_save, n_jobs, known_repeat = "CCCTAA")