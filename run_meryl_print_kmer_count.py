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

def run(path_meryldb, path_save, dir_tmp, dir_shell, dir_log, jobname, n_jobs, memory = 3000, hostname = list(), jobs_wait = list()):    
    cmd = f"""{PATH_MERYL} print threads={n_jobs} memory={memory} {path_meryldb} > {path_save}"""
        
    os.makedirs(os.path.dirname(path_save), exist_ok=True)
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
    
    list_meryl_count_target = [
        (
            "U10K_00016_hapmer",
            "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/PacBio/Hapmer_Diff/Hapmer_diff_k31/MATERNAL_Trimmed_F15_T5.noNbase.merged.k31.hapmer.meryl",
            "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/Binning_Compare/MATERNAL_Trimmed_F15_T5.noNbase.merged.k31.hapmer.kmer_count.txt"
        ),
        (
            "U10K_00020_hapmer",
            "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/PacBio/Hapmer_Diff/Hapmer_diff_k31/PATERNAL_Trimmed_F15_T5.noNbase.merged.k31.hapmer.meryl",
            "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/Binning_Compare/PATERNAL_Trimmed_F15_T5.noNbase.merged.k31.hapmer.kmer_count.txt"
        ),
        # (
        #     "U10K_00016_U10K_00020_Shared",
        #     "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/PacBio/Hapmer_Diff/Hapmer_diff_k31/PATERNAL_Trimmed_F15_T5.noNbase.merged.k31_and_MATERNAL_Trimmed_F15_T5.noNbase.merged.k31.meryl",
        #     "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/Binning_Compare/PATERNAL_Trimmed_F15_T5.noNbase.merged.k31_and_MATERNAL_Trimmed_F15_T5.noNbase.merged.k31.kmer_count.txt"
        # )
    ]
    
    curr_time = int(time.time())
    for name, path_meryldb, path_save in list_meryl_count_target:
        
        curr_time = int(time.time())
        run(
            path_meryldb, 
            path_save,
            f"{dir_workspace}/Meryl_Print_Count/tmp",
            f"{dir_workspace}/Meryl_Print_Count/sh",
            f"{dir_workspace}/Meryl_Print_Count/log",
            f"Meryl_Print_Count.{name}.{curr_time}",
            32,
            memory = 3000,
            hostname = [],
            jobs_wait = []
        )