#%%
import math, sys, os
from pathlib import Path

import pandas as pd
import numpy as np

sys.path.append(str(Path(__file__).parents[1]))
from fasta_parser import read_fasta_as_dict, save_specific_seqnames_from_fasta_dict


def main(path_hamming, path_fasta, hamming, col_origin_count, col_other_count, path_output):
    table_hamming = pd.read_csv(path_hamming, sep = '\t')
    table_hamming["tot_ham"] = table_hamming[col_origin_count] + table_hamming[col_other_count]
    table_hamming["hamming_error"] = 1 - table_hamming[col_origin_count]/table_hamming["tot_ham"]

    table_hamming["hamming_error"] = table_hamming["hamming_error"].fillna(0)

    dict_fasta = read_fasta_as_dict(path_fasta)

    list_seqnames_under_thresh = table_hamming[table_hamming["hamming_error"] < hamming]["Contig"].to_list()
    list_seqnames_save = list(filter(lambda val: val in list_seqnames_under_thresh, dict_fasta.keys()))

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

        # Flye Binned with0 (HiFi)
        "Flye_Binned_HiFi_with0":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20.with_0/Hap{hap}/assembly.fasta",

        # Flye Binned (HiFi)
        "Flye_Binned_HiFi":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20/Hap{hap}/assembly.fasta",

        # Flye Binned (ONT)
        "Flye_Binned_ONT":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/ONT_UL.Filter_L50kb_QV20_Q30Perc80/Hap{hap}/assembly.fasta",
        
        # # Flye Binned with 0 (ONT)
        # "Flye_Binned_ONT.with0":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/ONT_UL.Filter_L50kb_QV20_Q30Perc80.with_0/Hap{hap}/assembly.fasta",

        # # Unplaced Scaffold
        # "Unplaced_Scaffold" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/Chromosome_level_scaffolds.Assign_by_CHM13/CuratedDraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.hap{hap}.RagTag_with_Telomere_Assign.minimap2_latest.CHM13.unplaced.fasta"

    }
        
    dict_hap_to_parental = {
        1:'p',
        2:'m'
    }

    dict_hap_to_origin_col = {
        1:{
            "this":"PATERNAL_Trimmed_F15_T5.noNbase.merged.k21.hapmer",
            "other":"MATERNAL_Trimmed_F15_T5.noNbase.merged.k21.hapmer"
        },
        2:{
            "this":"MATERNAL_Trimmed_F15_T5.noNbase.merged.k21.hapmer",
            "other":"PATERNAL_Trimmed_F15_T5.noNbase.merged.k21.hapmer"
        }
    }

    hamming_thresh = 0.05
    
    for qv_filt in [70, 80, 90]:
        for hap in [1, 2]:
            for asm_name, path_asm_format in dict_path_assembly_format.items():
                path_asm = path_asm_format.format(hap = hap, triobin = dict_hap_to_parental[hap])

                path_asm_input = path_asm.replace(".fasta", f".qv{qv_filt}.fasta")

                assert os.path.exists(path_asm_input), path_asm_input

                dir_asm = os.path.dirname(path_asm)
                filename_asm = os.path.basename(path_asm)

                path_hamming = os.path.join(dir_asm, "Merqury_Hamming", f"{filename_asm.replace('.fasta', '')}.hapmers.count")
                assert os.path.exists(path_hamming), path_hamming

                path_output = os.path.join(
                    dir_asm,
                    filename_asm.replace(".fasta", f".qv{qv_filt}.hamming{hamming_thresh}.fasta")
                )

                main(path_hamming, path_asm_input, hamming_thresh, dict_hap_to_origin_col[hap]["this"], dict_hap_to_origin_col[hap]["other"], path_output)
                print(hap, asm_name)