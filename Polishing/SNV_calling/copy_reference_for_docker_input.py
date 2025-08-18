#%%
import os, subprocess

list_hap = [1,2]
list_triobins = ['m', 'p']

path_reference_format = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Reference_Fasta/Round1.SV_Polished_Fasta/PolishedASM.KOREF_S1.Trio_Specific.UpperCase.Polish_R1.SV.hap{hap}.fasta"

dir_output_format = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Read_Alignment/Align.Winnowmap2/Align_Polished_ASM.From_R1/Primary_Alignment/Docker_Input/ONT.Hap{hap}.triobin_{triobins}"


for hap in list_hap:
    for triobins in list_triobins:
        subprocess.run(
            f"cp {path_reference_format.format(hap = hap)} {dir_output_format.format(hap = hap, triobins = triobins)}", shell = True
        )
# %%
