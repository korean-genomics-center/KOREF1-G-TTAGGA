import pandas as pd
import os, time
import subprocess

PATH_MERYL = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/meryl-1.4.1/bin/meryl"


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

def run(kmer_size, list_reads, path_save_bin, dir_tmp, dir_shell, dir_log, jobname, n_jobs, memory = 3000, hostname = list(), jobs_wait = list()):    
    if len(list_reads) == 1:
        cmd = f"""{PATH_MERYL} k={kmer_size} count compress threads={n_jobs} memory={memory} output {path_save_bin} {list_reads[0]}"""
    
    else:
        list_tmp_counts = list()
        cmd = ""
        for ind, path_read in enumerate(list_reads, start = 1):
            path_save_count_of_read = f"{path_save_bin}.tmp.read{ind}"
            cmd += f"{PATH_MERYL} k={kmer_size} count compress threads={n_jobs} memory={memory} output {path_save_count_of_read} {path_read};\n"
            list_tmp_counts.append(path_save_count_of_read)
        cmd += f"{PATH_MERYL} union-sum threads={n_jobs} memory={memory} output {path_save_bin} {' '.join(list_tmp_counts)}"
        
    os.makedirs(os.path.dirname(path_save_bin), exist_ok=True)
    os.makedirs(dir_tmp, exist_ok=True)
    os.makedirs(dir_shell, exist_ok=True)
    os.makedirs(dir_log, exist_ok=True)
    
    write_shell_script_and_run(
        f"{dir_shell}/{jobname}.sh",
        get_context_for_qsub_submission(
            cmd,
            jobname,
            f"{dir_log}/{jobname}.stderr",
            f"{dir_log}/{jobname}.stdout",
            n_jobs,
            hostname=hostname,
            hold_jids=jobs_wait
        )
    )



if __name__ == "__main__":
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis"
    
    kmer_size = 30
    list_meryl_count_target = [
        (
            "U10K_00016",
            ["/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/MATERNAL/Trim_FQ/Merge/MATERNAL_Trimmed_F15_T5.noNbase.merged.R*.fq.gz"],
            f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/MATERNAL/Meryl_Binning/MATERNAL_Trimmed_F15_T5.noNbase.merged.R12.k{kmer_size}.compress.meryl"
        ),
        (
            "U10K_00020",
            ["/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/PATERNAL/Trim_FQ/Merge/PATERNAL_Trimmed_F15_T5.noNbase.merged.R*.fq.gz"],
            f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/PATERNAL/Meryl_Binning/PATERNAL_Trimmed_F15_T5.noNbase.merged.R12.k{kmer_size}.compress.meryl"
        ),
        # (
        #     "U10K_00016",
        #     ["/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/MATERNAL/Trim_FQ/Merge/MATERNAL_Trimmed_F15_T5.noNbase.merged.R1.fq.gz", "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/MATERNAL/Trim_FQ/Merge/MATERNAL_Trimmed_F15_T5.noNbase.merged.R2.fq.gz"],
        #     f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/MATERNAL/Meryl_Binning/MATERNAL_Trimmed_F15_T5.noNbase.merged.k{kmer_size}.compress.meryl"
        # ),
        # (
        #     "U10K_00020",
        #     ["/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/PATERNAL/Trim_FQ/Merge/PATERNAL_Trimmed_F15_T5.noNbase.merged.R1.fq.gz", "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/PATERNAL/Trim_FQ/Merge/PATERNAL_Trimmed_F15_T5.noNbase.merged.R2.fq.gz"],
        #     f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/PATERNAL/Meryl_Binning/PATERNAL_Trimmed_F15_T5.noNbase.merged.k{kmer_size}.compress.meryl"
        # ),
        # (
        #     "KOREF_HiFi",
        #     ["/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/PacBio/Merged/PacBio_HiFi_ccs_converted.filt.Filter_L1000_QV20.fastq.gz"],
        #     f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/PacBio/Meryl_Binning/PacBio_HiFi_ccs_converted.filt.Filter_L1000_QV20.k{kmer_size}.compress.meryl"
        # )
    ]
    
    curr_time = int(time.time())
    for name, list_read, path_save in list_meryl_count_target:
        
        curr_time = int(time.time())
        run(
            kmer_size,
            list_read,
            path_save,
            f"{dir_workspace}/Meryl_Binning/tmp",
            f"{dir_workspace}/Meryl_Binning/sh",
            f"{dir_workspace}/Meryl_Binning/log",
            f"Meryl_Binning.{name}.kmer{kmer_size}.compress.{curr_time}",
            32,
            memory = 800,
            hostname = [],
            jobs_wait = []
        )