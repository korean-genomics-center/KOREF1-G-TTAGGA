#%%
from pathlib import Path

import os, sys

sys.path.append(str(Path(__file__).parents[1] / "common"))
from fasta_parser import *
# %%

def run(path_asm1, path_asm2, path_save, list_chrnames_diploid, list_chrnames_haploid):
    dict_fasta_asm1 = read_fasta_as_dict(path_asm1)
    dict_fasta_asm2 = read_fasta_as_dict(path_asm2)

    dict_diploid = dict()
    for chrname in list_chrnames_diploid:
        dict_diploid[f"{chrname}_Hap1"] = dict_fasta_asm1[chrname]
        dict_diploid[f"{chrname}_Hap2"] = dict_fasta_asm2[chrname]
    
    for chrname in list_chrnames_haploid:
        in_hap1 = chrname in dict_fasta_asm1.keys()
        in_hap2 = chrname in dict_fasta_asm2.keys()

        assert not (in_hap1 == in_hap2), f"Not haploid: {chrname}"

        if in_hap1:
            dict_diploid[chrname] = dict_fasta_asm1[chrname]
        elif in_hap2:
            dict_diploid[chrname] = dict_fasta_asm2[chrname]
        else:
            raise Exception(chrname)
    
    print('\n'.join(list(dict_diploid.keys())))

    save_dict_fasta_into_file(dict_diploid, path_save)


#%%
if __name__ == "__main__":
    path_asm1 = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Final_Assembly/FinalASM.KOREF_S1.Trio_Specific.hap1.Reordered.fasta"
    path_asm2 = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Final_Assembly/FinalASM.KOREF_S1.Trio_Specific.hap2.Reordered.fasta"

    path_save = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Flagger_Evaluation/Diploid_Fasta/FinalASM.KOREF_S1.Trio_Specific.Diploid.Reordered.fasta"

    list_chrnames_diploid = list(map(lambda val: f"chr{val}", range(1, 23)))
    list_chrnames_haploid = ["chrX", "chrY"]

    run(path_asm1, path_asm2, path_save, list_chrnames_diploid, list_chrnames_haploid)
# %%
