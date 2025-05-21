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

def run_fastp_report(path_read1, path_read2, path_save_json, path_save_html, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):
    cmd_filtervariant = f"""{PATH_FASTP} --in1 {path_read1} --in2 {path_read2} --json {path_save_json} --html {path_save_html}"""
    
    for param_name, param_value in additional_params.items():
        cmd_filtervariant += f" --{param_name} {param_value}"
    
    os.makedirs(os.path.dirname(path_save_json), exist_ok=True)
    os.makedirs(os.path.dirname(path_save_html), exist_ok=True)
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



# if __name__ == "__main__":
#     dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis/Check_Trio_QC"
    
    
    
#     curr_time = int(time.time())
#     for sampleid in ["MATERNAL", "PATERNAL"]:
#         # for run in range(1, 9):
    
#         path_fasta_read1 = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/{sampleid}/Trim_FQ/Merge/{sampleid}_Trimmed_F15_T5.noNbase.merged.R1.fq.gz"
#         path_fasta_read2 = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/{sampleid}/Trim_FQ/Merge/{sampleid}_Trimmed_F15_T5.noNbase.merged.R2.fq.gz"
        
#         path_save_json = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/{sampleid}/Trim_FQ/Merge/Fastp_QC/{sampleid}_Trimmed_F15_T5.noNbase.merged.report.json"
#         path_save_html = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/{sampleid}/Trim_FQ/Merge/Fastp_QC/{sampleid}_Trimmed_F15_T5.noNbase.merged.report.html"
        
#         run_fastp_report(
#             path_fasta_read1,
#             path_fasta_read2,
#             path_save_json,
#             path_save_html,
#             f"{dir_workspace}/Check_QC/tmp",
#             f"{dir_workspace}/Check_QC/sh",
#             f"{dir_workspace}/Check_QC/log",
#             f"Report_Fastq.KOREF_Trio.{sampleid}.merged.{curr_time}",
#             16,
#             [],
#             []
#         )

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
        
        path_save_json = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/HiC/Deduplication/Fastp_Report/{filename_prefix}.dedup.report.json"
        path_save_html = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/HiC/Deduplication/Fastp_Report/{filename_prefix}.dedup.report.html"
        
        run_fastp_report(
            path_fasta_read1,
            path_fasta_read2,
            path_save_json,
            path_save_html,
            f"{dir_workspace}/Check_HiC_QC/tmp",
            f"{dir_workspace}/Check_HiC_QC/sh",
            f"{dir_workspace}/Check_HiC_QC/log",
            f"Report_Fastq.KOREF_HiC.{filename_prefix}.dedup.{curr_time}",
            16,
            [],
            []
        )