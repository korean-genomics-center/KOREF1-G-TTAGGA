#%%
import os, sys


sys.path.append("/BiO/Research/Project2/Endangered_Goral/Resources/Scripts/KOGIC_KOREF_AGGA/Gap_Filling/Manual_filling")
sys.path.append("/BiO/Research/Project2/Endangered_Goral/Resources/Scripts/KOGIC_KOREF_AGGA")
from check_contig_overlap_elongation import *
from fasta_parser import *

#%%
def read_fastq_as_dict(path_fastq):
    dict_fastq = dict()
    
    with open(path_fastq, 'r') as fr:
        while 1:
            header = fr.readline()
            if not header:break
            header = header.strip('\n').split()[0][1:]
            sequence = fr.readline().strip('\n')
            _ = fr.readline()
            qual = fr.readline().strip('\n')
            dict_fastq[header] = sequence
    return dict_fastq
#%%
dict_contig_to_seletion = {
    "h1tg000011l":{
        "select":"7bce8006-7b25-443d-a8e0-13ac5c238612",
        "revcomp":False
    },
    "h1tg000012l":{
        "select":"1cd42a03-8b01-4bc3-9215-3b97d6a713bf",
        "revcomp":False
    },
    "h1tg000026l":{
        "select":"c803a76f-bb59-45fb-8bb9-2b56f1bab222",
        "revcomp":True
    }
}

nround = 1

path_fastq_format = "/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/From_ONT/RawReads_Aligned_ONT/from_contig/{contig}_end.concise2.ont_readid.fastq"

path_save_fasta = "/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/From_ONT/Reads_for_Elongate/from_contig/{contig}_elongate_{nround}.fa"

# %%
def get_reverse_complement(sequence):
    dict_conv = {
        'A':'T',
        'T':'A',
        'G':'C',
        'C':'G'
    }
    
    comp_sequence = list(map(dict_conv.__getitem__, sequence))
    revcomp_sequence = ''.join(comp_sequence[::-1])
    
    return revcomp_sequence

for contig in dict_contig_to_seletion.keys():
    path_fastq = path_fastq_format.format(contig = contig)
    path_save = path_save_fasta.format(contig = contig, nround = nround)
    
    dict_fastq = read_fastq_as_dict(path_fastq)
    sequence_store = dict_fastq[dict_contig_to_seletion[contig]["select"]]
    
    if dict_contig_to_seletion[contig]["revcomp"]:
        sequence_store = get_reverse_complement(sequence_store)

    dict_store = {
        f"{contig}_elongate_{nround}":sequence_store
    }

    save_dict_fasta_into_file(dict_store, path_save)
# %%
