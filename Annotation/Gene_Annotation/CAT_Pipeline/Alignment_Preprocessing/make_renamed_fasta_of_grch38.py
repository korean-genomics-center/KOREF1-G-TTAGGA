#%%
import sys, os
from pathlib import Path

sys.path.append(str(Path(__file__).parents[4] / "common"))
from fasta_parser import *
# %%
path_grch38 = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/GRCh38/NCBI_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna"
path_save_renamed_primary = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Cactus_Alignment/GCF_000001405.40_GRCh38.p14_genomic.major_chromosome.id_conv.fa"


#%%
dict_fasta = read_fasta_as_dict(path_grch38)
# %%
list_chromosomes_major = list(filter(lambda val: val.startswith("NC"), dict_fasta.keys()))
# %%
dict_chrname_conv = dict()

for chrname in list_chromosomes_major:
    if "chromosome" not in chrname:
        continue
    chrname_suff = chrname.split(',')[0].split()[-1]
    
    chrname_conv = f"chr{chrname_suff}"
    
    dict_chrname_conv[chrname] = chrname_conv
    print(chrname_conv)
# %%
dict_renamed_fasta = dict()
for chrname, chrname_conv in dict_chrname_conv.items():
    dict_renamed_fasta[chrname_conv] = dict_fasta[chrname]
    
# %%
save_dict_fasta_into_file(dict_renamed_fasta, path_save_renamed_primary)
# %%
