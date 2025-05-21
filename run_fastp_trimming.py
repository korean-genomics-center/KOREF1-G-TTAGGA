import pandas as pd
import os, time
import subprocess

PATH_FASTP = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/Fastp-0.23.4/fastp"


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

def run_fastp_report(path_read1, path_read2, path_save_read1, path_save_read2, path_save_json, path_save_html, trim_front, trim_tail, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):
    cmd_filtervariant = f"""{PATH_FASTP} --in1 {path_read1} --in2 {path_read2} --out1 {path_save_read1} --out2 {path_save_read2} --json {path_save_json} --html {path_save_html} --n_base_limit 1 --detect_adapter_for_pe --cut_mean_quality 20 --average_qual 20 --length_required 15 --trim_front1 {trim_front} --trim_front2 {trim_front} --trim_tail1 {trim_tail} --trim_tail2 {trim_tail} --cut_front 1 --cut_tail 1 --cut_right 1 --thread {n_jobs}"""
    
    for param_name, param_value in additional_params.items():
        cmd_filtervariant += f" --{param_name} {param_value}"
    
    os.makedirs(os.path.dirname(path_save_read1), exist_ok=True)
    os.makedirs(os.path.dirname(path_save_read2), exist_ok=True)
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
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis"
    curr_time = int(time.time())
    for filename_prefix in [
"KOREF-S1_PBMC_Hi-C_HiSeqXten_20150518_DpnII-HinfI_L01",
"KOREF-S1_PBMC_Hi-C_HiSeqXten_20150518_MboI_L01",
"KOREF-S1_PBMC_Hi-C_NovaSeq_20190425_DpnII-HinfI_L01"
    ]:
    
        path_fasta_read1 = os.path.join("/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/HiC/Deduplication", f"{filename_prefix}_R1.dedup.fq.gz")
        path_fasta_read2 = os.path.join("/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/HiC/Deduplication", f"{filename_prefix}_R2.dedup.fq.gz")
        
        path_save_read1 = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/HiC/Trim_FQ/{filename_prefix}_R1.dedup.trim_F15_T0.fq.gz"
        path_save_read2 = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/HiC/Trim_FQ/{filename_prefix}_R2.dedup.trim_F15_T0.fq.gz"
        
        path_save_json = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/HiC/Trim_FQ/{filename_prefix}.dedup.trim_F15_T0.report.json"
        path_save_html = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/HiC/Trim_FQ/{filename_prefix}.dedup.trim_F15_T0.report.html"
        
        trim_front = 15
        trim_tail = 0
        
        curr_time = int(time.time())
        run_fastp_report(
            path_fasta_read1,
            path_fasta_read2,
            path_save_read1,
            path_save_read2,
            path_save_json,
            path_save_html,
            trim_front,
            trim_tail,
            f"{dir_workspace}/Trim_Fastq/tmp",
            f"{dir_workspace}/Trim_Fastq/sh",
            f"{dir_workspace}/Trim_Fastq/log",
            f"Trim_Fastq.KOREF_HiC.{filename_prefix}.dedup.Trim_F15_T0.{curr_time}",
            32,
            [],
            []
        )


if __name__ == "__main2__":
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis/Trim_Trio_FQ"
    
    
    curr_time = int(time.time())
    for sampleid in ["MATERNAL", "PATERNAL"]:
        for run in range(1, 9):
    
            path_fasta_read1 = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/KOREF-S1_Trio/{sampleid}/Blood/WGS/{sampleid}_L{run}_R1.fq.gz"
            path_fasta_read2 = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/KOREF-S1_Trio/{sampleid}/Blood/WGS/{sampleid}_L{run}_R2.fq.gz"
            
            path_save_read1 = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/{sampleid}/Trim_FQ/{sampleid}_L{run}_Trimmed_F0_T0.noNbase.R1.fq.gz"
            path_save_read2 = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/{sampleid}/Trim_FQ/{sampleid}_L{run}_Trimmed_F0_T0.noNbase.R2.fq.gz"
            
            path_save_json = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/{sampleid}/Trim_FQ/{sampleid}_L{run}_Trimmed_F0_T0.noNbase.report.json"
            path_save_html = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/{sampleid}/Trim_FQ/{sampleid}_L{run}_Trimmed_F0_T0.noNbase.report.html"
            
            trim_front = 0
            trim_tail = 0
            
            curr_time = int(time.time())
            run_fastp_report(
                path_fasta_read1,
                path_fasta_read2,
                path_save_read1,
                path_save_read2,
                path_save_json,
                path_save_html,
                trim_front,
                trim_tail,
                f"{dir_workspace}/Trim_Fastq/tmp",
                f"{dir_workspace}/Trim_Fastq/sh",
                f"{dir_workspace}/Trim_Fastq/log",
                f"Trim_Fastq.KOREF_Trio.{sampleid}.L{run}.{curr_time}",
                16,
                [],
                []
            )