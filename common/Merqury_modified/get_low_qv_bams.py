#%%
import math, os, subprocess
from joblib import Parallel, delayed

import pandas as pd
import numpy as np

PATH_SAMTOOLS = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/Samtools/samtools-1.17/samtools-1.17/bin/samtools"


def main(path_qv, path_bam, qv_thresh, output_prefix, n_jobs):
    table_qv = pd.read_csv(path_qv, names = ["seqname", "wrong", "total", "qv", "error"], sep = '\t')

    table_qv["chrname"] = table_qv["seqname"].apply(lambda val: val.split(':')[-2])
    table_qv["range_start"] = table_qv["seqname"].apply(lambda val: int(val.split(':')[-1].split('-')[0]))
    table_qv["range_end"] = table_qv["seqname"].apply(lambda val: int(val.split(':')[-1].split('-')[1]))

    table_lowqv = table_qv[table_qv["qv"] <= qv_thresh].reset_index(drop = True)

    with Parallel(n_jobs = n_jobs) as parallel:
        parallel(
            delayed(get_bam_of_region)(
                path_bam, output_prefix,
                table_lowqv.loc[ind, "chrname"],
                table_lowqv.loc[ind, "range_start"],
                table_lowqv.loc[ind, "range_end"]
            ) for ind in list(table_lowqv.index)
        )

def get_bam_of_region(path_bam, output_prefix, chrname, start, end):
    path_save = f"{output_prefix}.region.{chrname}.{start}_{end}.bam"
    cmd = f"{PATH_SAMTOOLS} view -h -b -@ 1 {path_bam} {chrname}:{start}-{end} > {path_save}"

    subprocess.run(cmd, shell = True)
    index_bam(path_save)

def index_bam(path_bam):
    cmd = f"{PATH_SAMTOOLS} index -@ 1 {path_bam}"
    subprocess.run(cmd, shell = True)


if __name__ == "__main__":
    path_bam_format = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Read_Alignment/Align.Winnowmap2/Align_Polished_ASM.From_R2/Primary_Alignment/Winnowmap_Align.Primary.Ref.KOREF_PolishedASM.R2_SNV.hap{hap}.QRY.PacBio_HiFi_ccs_converted.Filter_L1000_QV20.triobin_split.triobin_{triobin}0a.sorted.bam"
    path_qv_format = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Reference_Fasta/Round2.SNV_Polished_Fasta/Merqury/PolishedASM.KOREF_S1.Polish_R2.SNV.QV_per_window.k21/PolishedASM.KOREF_S1.Polish_R2.SNV.Hap{hap}.QV_per_window.k21.chunked.100000.qv"

    output_prefix_format = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Reference_Fasta/Round2.SNV_Polished_Fasta/Merqury/PolishedASM.KOREF_S1.Polish_R2.SNV.QV_per_window.k21/Align_on_LowQV_Region/Hap{hap}.PacBio_HiFi_ccs_converted.Filter_L1000_QV20.triobin_{triobin}0a/Winnowmap_Align.Primary.Ref.KOREF_PolishedASM.R2_SNV.hap{hap}.QRY.PacBio_HiFi_ccs_converted.Filter_L1000_QV20.triobin_split.triobin_{triobin}0a.sorted"

    dict_hap_to_triobin = {
        1:'p',
        2:'m'
    }

    qv_thresh = 40
    n_jobs = 10

    for hap in [1, 2]:
        path_bam = path_bam_format.format(
            hap = hap, triobin = dict_hap_to_triobin[hap]
        )
        path_qv = path_qv_format.format(hap = hap)
        output_prefix = output_prefix_format.format(hap = hap, triobin = dict_hap_to_triobin[hap])

        assert os.path.exists(path_bam)
        assert os.path.exists(path_qv)
        
        os.makedirs(os.path.dirname(output_prefix), exist_ok=True)

        main(path_qv, path_bam, qv_thresh, output_prefix, n_jobs)