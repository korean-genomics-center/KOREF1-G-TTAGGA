#%%
import os
from copy import deepcopy

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

def run(list_path_salmon_qf, tpm_cutoff, path_output_prefix, draw_hist = True):
    os.makedirs(os.path.dirname(path_output_prefix), exist_ok=True)

    list_transcripts_unexpressed = list()
    for ind_file, path_salmon_qf in enumerate(list_path_salmon_qf, start = 1):
        table_salmon_qf = pd.read_csv(path_salmon_qf, sep = '\t')

        list_tpms = table_salmon_qf["TPM"].to_list()

        if draw_hist:
            draw_histogram(list_tpms, f"{path_output_prefix}.tpm_hist.file_{ind_file}.png", color = "gray", bins = 100)
            draw_histogram(list_tpms, f"{path_output_prefix}.tpm_hist.file_{ind_file}.bin_to10.png", color = "gray", bins = np.linspace(0, 10, 101))
            draw_histogram(list_tpms, f"{path_output_prefix}.tpm_hist.file_{ind_file}.logscale.png", color = "gray", bins = 100, yscale = "log")
            draw_histogram(list_tpms, f"{path_output_prefix}.tpm_hist.file_{ind_file}.bin_to10.logscale.png", color = "gray", bins = np.linspace(0, 10, 101), yscale = "log")
            draw_histogram(list_tpms, f"{path_output_prefix}.tpm_hist.file_{ind_file}.bin_to5.png", color = "gray", bins = np.linspace(0, 5, 101))
            draw_histogram(list_tpms, f"{path_output_prefix}.tpm_hist.file_{ind_file}.bin_to5.logscale.png", color = "gray", bins = np.linspace(0, 5, 101), yscale = "log")
            draw_histogram(list_tpms, f"{path_output_prefix}.tpm_hist.file_{ind_file}.bin_to2.png", color = "gray", bins = np.linspace(0, 2, 101))
            draw_histogram(list_tpms, f"{path_output_prefix}.tpm_hist.file_{ind_file}.bin_to2.logscale.png", color = "gray", bins = np.linspace(0, 2, 101), yscale = "log")
        
        list_transcripts_under_threshold_thisfile = table_salmon_qf[table_salmon_qf["TPM"] < tpm_cutoff]["Name"].to_list()

        if len(list_transcripts_unexpressed) == 0:
            list_transcripts_unexpressed = deepcopy(list_transcripts_under_threshold_thisfile)
        else:
            list_transcripts_unexpressed = list(set(list_transcripts_unexpressed) & set(list_transcripts_under_threshold_thisfile))
        print(len(list_transcripts_unexpressed))

    path_save_transcript_ids = f"{path_output_prefix}.TPM_cutoff_{tpm_cutoff}.transcripts_notexpressed.txt"
    with open(path_save_transcript_ids, 'w') as fw:
        for transcript in list_transcripts_unexpressed:
            fw.write(transcript+'\n')


def draw_histogram(list_tpm, path_save, yscale = "linear", **kwargs):
    plt.rcParams["font.size"] = 14
    fig, ax = plt.subplots(1, 1, figsize = (4, 3))
    ax.hist(list_tpm, **kwargs)
    ax.set_xlabel("TPM")
    ax.set_ylabel("# Transcripts")
    ax.set_yscale(yscale)
    plt.savefig(path_save, facecolor = 'w', bbox_inches = "tight")
    plt.close()



if __name__ == "__main__":

    for hap in [1, 2]:
        for tpm_cutoff in [0.5, 1]:
            list_path_salmon_qf = [
                f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/CDS_Fasta/Hap{hap}.Blood_Cellline.CGP_PB/Salmon/Quant_Salmon/RNASeq_Blood_NovaSeq_20190403/quant.sf",
                f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/CDS_Fasta/Hap{hap}.Blood_Cellline.CGP_PB/Salmon/Quant_Salmon/RNASeq_CellLine_HiSeq_20180417/quant.sf"
            ]
            path_output_prefix = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/CDS_Fasta/Hap{hap}.Blood_Cellline.CGP_PB/Salmon/Quant_Salmon/Filter_by_TPM/Transcript_Filtering.Hap{hap}"

            run(list_path_salmon_qf, tpm_cutoff, path_output_prefix, draw_hist=tpm_cutoff==1)