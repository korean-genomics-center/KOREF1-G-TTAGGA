#%%
import math, sys, os
from pathlib import Path

import pandas as pd
import numpy as np

sys.path.append(str(Path(__file__).parents[1]))
from fasta_parser import read_fasta_as_dict, save_specific_seqnames_from_fasta_dict


def main(path_qv, path_fasta, qv, path_output):
    table_qv = pd.read_csv(path_qv, names = ["seqname", "wrong", "total", "qv", "error"], sep = '\t')

    dict_fasta = read_fasta_as_dict(path_fasta)

    list_seqnames_over_qv = table_qv[table_qv["qv"] > qv]["seqname"].to_list()
    list_seqnames_save = list(filter(lambda val: val in list_seqnames_over_qv, dict_fasta.keys()))
    print(len(list_seqnames_save), len(dict_fasta))
    save_specific_seqnames_from_fasta_dict(dict_fasta, list_seqnames_save, path_output)


if __name__ == "__main__":
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis"
    
    dict_path_assembly_format = {
        # Hifiasm Binned (HiFi)
        "Hifiasm_Binned_HiFi":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20/Hap{hap}/KOREF_S1.Binned_Assembly.hifi_L1000_QV20.triobin_{triobin}.haploid.l0.asm.bp.p_ctg.fasta",

        #Hifiasm Binned (Alt) (HiFi)
        "Hifiasm_Binned_ALT_HiFi":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20/Hap{hap}/KOREF_S1.Binned_Assembly.hifi_L1000_QV20.triobin_{triobin}.haploid.l0.asm.bp.a_ctg.fasta",

        # Hifiasm Binned with0 (HiFi)
        "Hifiasm_Binned_HiFi_with0":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20.with_0/Hap{hap}/KOREF_S1.Binned_Assembly.hifi_L1000_QV20.triobin_{triobin}_0.haploid.l0.asm.bp.p_ctg.fasta",
        
        #Hifiasm Binned with0 (Alt) (HiFi)
        "Hifiasm_Binned_ALT_HiFi_with0":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20.with_0/Hap{hap}/KOREF_S1.Binned_Assembly.hifi_L1000_QV20.triobin_{triobin}_0.haploid.l0.asm.bp.a_ctg.fasta",

        # Flye Binned (ONT)
        "Flye_Binned_ONT":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/ONT_UL.Filter_L50kb_QV20_Q30Perc80/Hap{hap}/assembly.fasta",
        
        # # Flye Binned with 0 (ONT)
        # "Flye_Binned_ONT.with0":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/ONT_UL.Filter_L50kb_QV20_Q30Perc80.with_0/Hap{hap}/assembly.fasta",

        # Flye Binned with0 (HiFi)
        "Flye_Binned_HiFi_with0":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20.with_0/Hap{hap}/assembly.fasta",

        # Flye Binned (HiFi)
        "Flye_Binned_HiFi":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20/Hap{hap}/assembly.fasta"
    }
        
    dict_hap_to_parental = {
        1:'p',
        2:'m'
    }

    qv_thresh = 80

    for hap in [1, 2]:
        for asm_name, path_asm_format in dict_path_assembly_format.items():
            path_asm = path_asm_format.format(hap = hap, triobin = dict_hap_to_parental[hap])

            assert os.path.exists(path_asm), path_asm

            dir_asm = os.path.dirname(path_asm)
            filename_asm = os.path.basename(path_asm)

            path_qv = os.path.join(dir_asm, "Merqury_QV", f"{filename_asm.replace('.fasta', '')}.k31.perseq.qv")

            path_output = os.path.join(
                dir_asm,
                filename_asm.replace(".fasta", f".qv{qv_thresh}.fasta")
            )

            main(path_qv, path_asm, qv_thresh, path_output)
            print(hap, asm_name)