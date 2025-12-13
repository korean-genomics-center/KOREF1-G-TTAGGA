#%%
import os, sys


sys.path.append("/BiO/Research/Project2/Endangered_Goral/Resources/Scripts/KOGIC_KOREF_AGGA/Gap_Filling/Manual_filling")
sys.path.append("/BiO/Research/Project2/Endangered_Goral/Resources/Scripts/KOGIC_KOREF_AGGA")
from check_contig_overlap_elongation import *
from fasta_parser import *


#%%
path_paf = "/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Assembly_Elongation/from_round1/Minimap2_Align/Align_Minimap2.QRY_LocalASM_from_h1tg000012l.round2.maphifiPcw5s100.paf"

table_paf = read_paf(path_paf)
# %%
table_paf_filt = filter_paf_of_end_aligned_query(table_paf, min_idy = 0.99, min_align = 1000, max_hang = 10000, min_align2hang = 0.99)
#%%

contig_select = "ptg000001l"
revcomp = False
nround = 2

path_fasta = "/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Assembly_Elongation/from_round1/h1tg000012l_elongate_1_end.Hifiasm.HiFi/h1tg000012l_elongate_1_end.hifi.haploid.l0.asm.bp.p_ctg.fasta"

path_save_fasta = f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Assembly_Elongation/from_round1/Final_contig_elongate/h1tg000012l_elongate_{nround}.fa"

dict_fasta = read_fasta_as_dict(path_fasta)
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

if revcomp:
    sequence_store = get_reverse_complement(
        dict_fasta[contig_select]
    )
else:
    sequence_store = dict_fasta[contig_select]
# %%
dict_store = {
    f"h1tg000012l_elongate_{nround}":sequence_store
}

save_dict_fasta_into_file(dict_store, path_save_fasta)
# %%
