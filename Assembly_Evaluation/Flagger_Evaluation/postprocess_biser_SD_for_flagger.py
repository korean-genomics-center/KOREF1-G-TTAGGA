#%%
import os,sys,math

import pandas as pd
import numpy as np


def run(path_sd, path_out, error_thresh):
    table_sd = pd.read_csv(path_sd, sep = '\t', names = list(range(14)), usecols = list(range(12)))
    table_sd = table_sd.rename(
        columns = dict(zip(
            list(range(12)),
            ["chr1", "start1", "end1", "chr2", "start2", "end2", "ref", "score", "strand1", "strand2", "maxlen", "alnlen", "cigar", "opt"]
        ))
    )

    table_sd_filt = table_sd[table_sd["score"] < error_thresh]
    
    table_sd_filt_interchr = table_sd_filt[table_sd_filt["chr1"] != table_sd_filt["chr2"]]
    table_sd_filt_intrachr = table_sd_filt[table_sd_filt["chr1"] == table_sd_filt["chr2"]]

    table_sd_filt_intrachr_filt = table_sd_filt_intrachr[~np.logical_and(
        table_sd_filt_intrachr["start1"] == table_sd_filt_intrachr["start2"],
        table_sd_filt_intrachr["end1"] == table_sd_filt_intrachr["end2"]
    )]

    table_sd_out = pd.concat([table_sd_filt_interchr, table_sd_filt_intrachr_filt])

    table_sd_out_bed1 = table_sd_out[["chr1", "start1", "end1"]].rename(columns = {"chr1":"chr", "start1":"start", "end1":"end"})
    table_sd_out_bed2 = table_sd_out[["chr2", "start2", "end2"]].rename(columns = {"chr2":"chr", "start2":"start", "end2":"end"})

    table_sd_out_bed = pd.concat([table_sd_out_bed1, table_sd_out_bed2]).drop_duplicates(subset = ["chr", "start", "end"])
    table_sd_out_bed.to_csv(path_out, sep = '\t', header = False, index = False)


#%%
if __name__ == "__main__":
    error_thresh = 10

    for hap in [1, 2]:
        path_sd = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Repeat_Annotation/Segmental_Duplication/Biser/Hap{hap}/FinalASM.KOREF_S1.Trio_Specific.hap{hap}.Reordered.masked.SD"
        path_out = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Flagger_Evaluation/Annotations/SD.hap{hap}.err_l{error_thresh}.bed"

        run(path_sd, path_out, error_thresh)