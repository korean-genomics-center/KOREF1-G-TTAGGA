import pandas as pd
import os, time
import subprocess

PATH_SCRIPT = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Resources/Scripts/Read_Filtering/merge_fqgz_files.py"


def get_context_for_qsub_submission(cmd, jobname, stderr_path, stdout_path, nthreads, hostname = list(), hold_jids = list()):
    #Header
    script_header = f"#!bin/sh\n"
    script_header += f"#$ -N {jobname}\n"
    script_header += f"#$ -S /bin/bash\n"
    script_header += f"#$ -V\n"
    if len(hold_jids) > 0:
        script_header += f"#$ -hold_jid {','.join(hold_jids)}\n"
    if len(hostname) > 0:
        script_header += f"#$ -l hostname={'|'.join(hostname)}\n"
    script_header += f"#$ -e {stderr_path}\n"
    script_header += f"#$ -o {stdout_path}\n"
    script_header += f"#$ -pe smp {nthreads}\n"
    script_header +="\n"
    script_header+=f"echo \"=================================================\"\n"
    script_header+=f"echo \"Started on	   : $(date)\"\n"
    script_header+=f"echo \"Running on node  : $(hostname)\"\n"
    script_header+=f"echo \"Current job ID   : $JOB_ID\"\n"
    script_header+=f"echo \"Current job name : $JOB_NAME\"\n"
    script_header+=f"echo \"=================================================\"\n\n"

    script_footer="\n\n"
    script_footer+=f"echo \"=================================================\"\n"
    script_footer+=f"echo \"Finished on	   : $(date)\"\n"
    script_footer+=f"echo \"=================================================\"\n"
    script_footer+=f"echo \"Done.\"\n"

    script = script_header + cmd + script_footer
    return script

def write_shell_script_and_run(path_shell, shell_context):
    with open(path_shell, 'w') as fw:
        fw.write(shell_context)
    subprocess.run(f"qsub {path_shell}", shell = True)

def run_merge(list_files, path_save, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):
    cmd_filtervariant = f"""python {PATH_SCRIPT} --files {' '.join(list_files)} --out {path_save}"""
    
    for param_name, param_value in additional_params.items():
        cmd_filtervariant += f" --{param_name} {param_value}"
    
    os.makedirs(os.path.dirname(path_save), exist_ok=True)
    os.makedirs(dir_tmp, exist_ok=True)
    os.makedirs(dir_shell, exist_ok=True)
    os.makedirs(dir_log, exist_ok=True)
    
    write_shell_script_and_run(
        f"{dir_shell}/{jobname}.sh",
        get_context_for_qsub_submission(
            cmd_filtervariant,
            jobname,
            f"{dir_log}/{jobname}.stderr",
            f"{dir_log}/{jobname}.stdout",
            n_jobs,
            hostname=hostname,
            hold_jids=jobs_wait
        )
    )

if __name__ == "__main__":
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis/Merge_FQ"
    
    
    curr_time = int(time.time())

    
    for pair in [('p', 'a', '0')]:
        list_fasta_read = list(map(lambda triobin: f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/PacBio/Parental_Split/PacBio_HiFi_ccs_converted.filt.Filter_L1000_QV20.merged.triobin_split.triobin_{triobin}.fq.gz", pair))
    
        path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/PacBio/Parental_Split/PacBio_HiFi_ccs_converted.filt.Filter_L1000_QV20.merged.triobin_split.triobin_{''.join(pair)}.fq.gz"
        
        curr_time = int(time.time())
        run_merge(
            list_fasta_read,
            path_save,
            f"{dir_workspace}/Merge_Fastq/tmp",
            f"{dir_workspace}/Merge_Fastq/sh",
            f"{dir_workspace}/Merge_Fastq/log",
            f"Merge_Fastq.KOREF-S1_HiFi.{''.join(pair)}.{curr_time}",
            1,
            [],
            []
        )

if __name__ == "__main3__":
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis/Merge_FQ"
    
    
    curr_time = int(time.time())
    for sampleid in ["MATERNAL", "PATERNAL"]:
        for pair  in [1, 2]:
            list_fasta_read = list(map(lambda run: f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/{sampleid}/Trim_FQ/{sampleid}_L{run}_Trimmed_F15_T5.noNbase.R{pair}.fq.gz", range(1, 9)))
        
            path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/{sampleid}/Trim_FQ/Merge/{sampleid}_Trimmed_F15_T5.noNbase.merged.R{pair}.fq.gz"
            
            curr_time = int(time.time())
            run_merge(
                list_fasta_read,
                path_save,
                f"{dir_workspace}/Merge_Fastq/tmp",
                f"{dir_workspace}/Merge_Fastq/sh",
                f"{dir_workspace}/Merge_Fastq/log",
                f"Merge_Fastq.KOREF_Trio.{sampleid}.R{pair}.{curr_time}",
                1,
                [],
                [
                    "Trim_Fastq.KOREF_Trio.MATERNAL.L1.1745692871",
                    "Trim_Fastq.KOREF_Trio.MATERNAL.L2.1745692871",
                    "Trim_Fastq.KOREF_Trio.MATERNAL.L3.1745692871",
                    "Trim_Fastq.KOREF_Trio.MATERNAL.L4.1745692871",
                    "Trim_Fastq.KOREF_Trio.MATERNAL.L5.1745692871",
                    "Trim_Fastq.KOREF_Trio.MATERNAL.L6.1745692871",
                    "Trim_Fastq.KOREF_Trio.MATERNAL.L7.1745692871",
                    "Trim_Fastq.KOREF_Trio.MATERNAL.L8.1745692871",
                    "Trim_Fastq.KOREF_Trio.PATERNAL.L1.1745692871",
                    "Trim_Fastq.KOREF_Trio.PATERNAL.L2.1745692871",
                    "Trim_Fastq.KOREF_Trio.PATERNAL.L3.1745692871",
                    "Trim_Fastq.KOREF_Trio.PATERNAL.L4.1745692871",
                    "Trim_Fastq.KOREF_Trio.PATERNAL.L5.1745692871",
                    "Trim_Fastq.KOREF_Trio.PATERNAL.L6.1745692872",
                    "Trim_Fastq.KOREF_Trio.PATERNAL.L7.1745692872",
                    "Trim_Fastq.KOREF_Trio.PATERNAL.L8.1745692872"
                ]
            )
            
if __name__ == "__main3__":
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis/Merge_FQ"
    
    list_files = [
        "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/PacBio/Adapter_Trimmed/KOREF-S1_Blood_WGS_Sequel-II_20190614_FC00/KOREF-S1_Blood_WGS_Sequel-II_20190614_FC00.ccs.converted.filt.Filter_L1000_QV20.fastq.gz",
        "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/PacBio/Adapter_Trimmed/KOREF-S1_Blood_WGS_Sequel-II_20190709_FC00/KOREF-S1_Blood_WGS_Sequel-II_20190709_FC00.ccs.converted.filt.Filter_L1000_QV20.fastq.gz",
        "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/PacBio/Adapter_Trimmed/KOREF-S1_Blood_WGS_Sequel-II_20190822_FC00/KOREF-S1_Blood_WGS_Sequel-II_20190822_FC00.ccs.converted.filt.Filter_L1000_QV20.fastq.gz",
        "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/PacBio/Adapter_Trimmed/KOREF-S1_Blood_WGS_Sequel-II_20190918_FC00/KOREF-S1_Blood_WGS_Sequel-II_20190918_FC00.ccs.converted.filt.Filter_L1000_QV20.fastq.gz",
        "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/PacBio/Adapter_Trimmed/KOREF-S1_Blood_WGS_Sequel-II_20190923_FC00/KOREF-S1_Blood_WGS_Sequel-II_20190923_FC00.ccs.converted.filt.Filter_L1000_QV20.fastq.gz",
        "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/PacBio/Adapter_Trimmed/KOREF-S1_Blood_WGS_Sequel-II_20191210_FC00/KOREF-S1_Blood_WGS_Sequel-II_20191210_FC00.ccs.converted.filt.Filter_L1000_QV20.fastq.gz",
        "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/PacBio/Adapter_Trimmed/TN2409D0573_CELL1/TN2409D0573_CELL1.hifi_reads.bc2073.converted.filt.Filter_L1000_QV20.fastq.gz",
        "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/PacBio/Adapter_Trimmed/TN2409D0573_CELL2/TN2409D0573_CELL2.hifi_reads.bc2073.converted.filt.Filter_L1000_QV20.fastq.gz"
    ]
    run_merge(
        list_files,
        "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/PacBio/Merged/PacBio_HiFi_ccs_converted.filt.Filter_L1000_QV20.fastq.gz",
        f"{dir_workspace}/Merge_Fastq/tmp",
        f"{dir_workspace}/Merge_Fastq/sh",
        f"{dir_workspace}/Merge_Fastq/log",
        f"Merge_Fastq.KOREF_HiFi.{int(time.time())}",
        1,
        [],
        []
    )
    
if __name__ == "__main2__":
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis/Merge_FQ"
    
    list_seqname = [
        "20240517_1158_P2S-01506-B_PAW01562_215af4e8",
        "20240629_1736_P2S-01506-A_PAW01687_bf35a486",
        "20240629_1736_P2S-01506-B_PAU99972_b4ed6da3",
        "20240811_1107_P2S-01506-A_PAW77551_0cd923fa",
        "20240816_0603_P2S-01501-B_PAW77472_d96fb6c9",
        "20240817_0645_P2S-01501-B_PAW77472_c8fb51b0",
        "20240818_1702_P2S-01506-B_PAW88380_5c7050f8",
        "20240920_1805_P2S-01506-A_PAY14781_7526b2f4",
        "20240923_1616_P2S-01506-A_PAY08433_74256187",
        "20240923_1641_P2S-01506-B_PAY03655_fb6f2a12",
        "20240924_2038_P2S-01506-A_PAY08433_7a23907e",
        "20240924_2038_P2S-01506-B_PAY03655_f02a3d27",
        "20240926_1205_P2S-01506-A_PAY09654_90f482e1",
        "20240926_1205_P2S-01506-B_PAY09817_4744eb86",
        "20240930_1539_P2S-01506-B_PAY10807_076f0b43",
        "20241002_1732_P2S-01506-A_PAY10807_3ab61abe",
        "20241011_1754_P2S-01506-A_PAW99288_c001806f",
        "20241012_1341_P2S-01506-A_PAW99288_c255cb62",
        "20241012_1341_P2S-01506-B_PAY01199_fef5a120",
        "20241021_1517_P2S-01506-A_PAW99288_e821305c",
        "20241021_1517_P2S-01506-B_PAY01199_b712d259",
        "20241022_1541_P2S-01506-A_PAW99013_6c642a53",
        "20241022_1541_P2S-01506-B_PAY07547_1caa91a9",
        "20241031_1212_P2S-01506-A_PAY83904_baf56bf4",
        "20241031_1220_P2S-01506-B_PAY74072_775e6f78",
        "20241031_1440_P2S-01503-A_PAY74208_5cb0c269",
        "20241031_1440_P2S-01503-B_PAY74099_bed3405c",
        "20241101_1326_P2S-01506-A_PAY83904_5bb2fcbf",
        "20241101_1326_P2S-01506-B_PAY74072_95afdfdf"
    ]
    
    list_files = list(map(lambda seqname: f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/ONT/UL_Filtered/{seqname}.Adapt_Filt.Edge_Trim_Front_1500_Tail_1500.Filter_L100kb_QV15.fastq.gz", list_seqname))
    
    run_merge(
        list_files,
        "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/ONT/Merged/ONT_UL.Adapt_Filt.Edge_Trim_Front_1500_Tail_1500.Filter_L100kb_QV15.merged.fastq.gz",
        f"{dir_workspace}/Merge_Fastq/tmp",
        f"{dir_workspace}/Merge_Fastq/sh",
        f"{dir_workspace}/Merge_Fastq/log",
        f"Merge_Fastq.KOREF_ONT.{int(time.time())}",
        1,
        [],
        []
    )