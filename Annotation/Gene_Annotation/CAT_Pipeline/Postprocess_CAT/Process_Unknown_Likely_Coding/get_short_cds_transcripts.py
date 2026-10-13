#%%
import os,sys
from copy import deepcopy

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from pathlib import Path

sys.path.append(str(Path(__file__).parents[4]))
from fasta_parser import read_fasta_as_dict


def run(path_cds_fasta, length_cutoff, path_output_prefix, draw_hist = True):
    os.makedirs(os.path.dirname(path_output_prefix), exist_ok=True)

    dict_cds_fasta = read_fasta_as_dict(path_cds_fasta)
    dict_cds_length = {cdsname:len(dict_cds_fasta[cdsname]) for cdsname in dict_cds_fasta.keys()}

    list_transcripts_short = list(filter(lambda name: dict_cds_length[name] < length_cutoff, dict_cds_length.keys()))
    
    if draw_hist:
        list_lengths = list(dict_cds_length.values())
        draw_histogram(list_lengths, f"{path_output_prefix}.cdslen_hist.png", color = "gray", bins = 100)
        draw_histogram(list_lengths, f"{path_output_prefix}.cdslen_hist.logscale.png", color = "gray", bins = 100, yscale = "log")
        draw_histogram(list_lengths, f"{path_output_prefix}.cdslen_hist.bin_to1k.png", color = "gray", bins = np.linspace(0, 1000, 101))
        draw_histogram(list_lengths, f"{path_output_prefix}.cdslen_hist.bin_to1k.logscale.png", color = "gray", bins = np.linspace(0, 1000, 101), yscale = "log")

    path_save_transcript_ids = f"{path_output_prefix}.CDSLength_cutoff_{length_cutoff}.tooshort.txt"
    with open(path_save_transcript_ids, 'w') as fw:
        for transcript in list_transcripts_short:
            fw.write(transcript+'\n')

def draw_histogram(list_length, path_save, yscale = "linear", **kwargs):
    plt.rcParams["font.size"] = 14
    fig, ax = plt.subplots(1, 1, figsize = (4, 3))
    ax.hist(list_length, **kwargs)
    ax.set_xlabel("CDS Length (bp)")
    ax.set_ylabel("# Transcripts")
    ax.set_yscale(yscale)
    plt.savefig(path_save, facecolor = 'w', bbox_inches = "tight")
    plt.close()


#%%
if __name__ == "__main__":
    length_cutoff = 100

    for hap in [1, 2]:
        for gff_type in ["whole", "known", "unknown"]:
            path_cds_fasta = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/CDS_Fasta/Hap{hap}.Blood_Cellline.CGP_PB/KOREF1_ASM_Hap{hap}.{gff_type}_gff.transcript_CDS.fa"
            path_output_prefix = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/CDS_Fasta/Hap{hap}.Blood_Cellline.CGP_PB/CDS_Length/GFF_{gff_type}/Transcript_Filtering.Hap{hap}.GFF_{gff_type}"

            run(path_cds_fasta, length_cutoff, path_output_prefix, draw_hist = True)