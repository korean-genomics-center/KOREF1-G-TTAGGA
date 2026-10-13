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

    list_runs = ["201904", "20210426"]
    
    for readpair in [1, 2]:
        list_fasta_read = list(map(lambda run: f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/IlluminaWGS/RemoveAdapter/KOREF1_WGS_NovaSeq_{run}_R{readpair}.trimPolyG.Trim_F15_T5.fq.gz", list_runs))
    
        path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/IlluminaWGS/Merged_Read/KOREF1_WGS_NovaSeq_ReadMerged_R{readpair}.trimPolyG.Trim_F15_T5.fq.gz"
        
        curr_time = int(time.time())
        run_merge(
            list_fasta_read,
            path_save,
            f"{dir_workspace}/Merge_Fastq/tmp",
            f"{dir_workspace}/Merge_Fastq/sh",
            f"{dir_workspace}/Merge_Fastq/log",
            f"Merge_Fastq.WGS_NovaSeq.R{readpair}.{curr_time}",
            1,
            [],
            []
        )