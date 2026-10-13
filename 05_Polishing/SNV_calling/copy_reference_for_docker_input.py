#%%
import os, subprocess

list_hap = [1,2]
#%%
path_reference_format = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Reference_Fasta/Round8.SV_Polished_Fasta/PolishedASM.KOREF_S1.Trio_Specific.UpperCase.Polish_R8.SV.hap{hap}.fasta"

dir_output_format = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Read_Alignment/Align.Winnowmap2/Align_PolishedASM.Polish_R8.SV/Primary_Alignment/Docker_Input/{tgs}.Hap{hap}.triobin_{triobins}"

dict_hap_to_parental = {
    1:"p",
    2:"m" 
}

dict_hap_exchange = {
    'p' : 'm',
    'm' : 'p'
}
#%%
for tgs in ["ONT", "HiFi"]:
    for hap in list_hap:
        if tgs == "HiFi":
            list_triobins = ["pm0a", f"{dict_hap_to_parental[hap]}0a", dict_hap_to_parental[hap], dict_hap_exchange[dict_hap_to_parental[hap]]]
        else:
            list_triobins = [f"{dict_hap_to_parental[hap]}0a", dict_hap_to_parental[hap], dict_hap_exchange[dict_hap_to_parental[hap]]]
        for triobins in list_triobins:
            dir_out = dir_output_format.format(tgs = tgs, hap = hap, triobins = triobins)
            os.makedirs(dir_out, exist_ok = True)
            subprocess.run(
                f"cp {path_reference_format.format(hap = hap)} {dir_out}", shell = True
            )
# %%
if __name__ == "__tmp__":
    path_bam_format = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Read_Alignment/Align.Winnowmap2/Align_PolishedASM.Polish_R1.SV/Primary_Alignment/Winnowmap_Align.Ref.KOREF_PolishedASM.Polish_R1.SV.hap{hap}.QRY.PacBio_HiFi_ccs_converted.Filter_L1000_QV20.triobin_split.triobin_{triobins}.sorted.bam"
    path_bai_format = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Read_Alignment/Align.Winnowmap2/Align_PolishedASM.Polish_R1.SV/Primary_Alignment/Winnowmap_Align.Ref.KOREF_PolishedASM.Polish_R1.SV.hap{hap}.QRY.PacBio_HiFi_ccs_converted.Filter_L1000_QV20.triobin_split.triobin_{triobins}.sorted.bam.bai"

    for hap in list_hap:
        for triobins in ["pm0a", f"{dict_hap_to_parental[hap]}0a",  dict_hap_to_parental[hap], dict_hap_exchange[dict_hap_to_parental[hap]]]:
            dir_out = dir_output_format.format(hap = hap, triobins = triobins, tgs = "HiFi")
            subprocess.run(
                f"mv {path_bam_format.format(hap = hap, triobins = triobins)} {dir_out}", shell = True
            )
            subprocess.run(
                f"mv {path_bai_format.format(hap = hap, triobins = triobins)} {dir_out}", shell = True
            )
# %%
