import pandas as pd
import os, time
import subprocess

PATH_MERYL_FOR_WINNNOWMAP = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/winnowmap/winnowmap.v2.03/Winnowmap/bin/meryl"


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

def run(kmer_size, path_fasta, path_save_prefix, dir_tmp, dir_shell, dir_log, jobname, n_jobs, memory = 3000, hostname = list(), jobs_wait = list()):    
    path_save_meryldb = f"{path_save_prefix}_MerylDB"
    path_save_merylprint = f"{path_save_prefix}.repetitive_k{kmer_size}.txt"
    cmd = f"""{PATH_MERYL_FOR_WINNNOWMAP} count k={kmer_size} threads={n_jobs} memory={memory} output {path_save_meryldb} {path_fasta}; {PATH_MERYL_FOR_WINNNOWMAP} print greater-than distinct=0.9998 {path_save_meryldb} > {path_save_merylprint}"""
        
    os.makedirs(os.path.dirname(path_save_prefix), exist_ok=True)
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
    
    kmer_size = 15
    
    path_fasta = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Flagger_Evaluation/Diploid_Fasta/FinalASM.KOREF_S1.Trio_Specific.Diploid.Reordered.fasta"
    path_save_prefix = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Flagger_Evaluation/Diploid_Fasta/MerylDB/KOREF_FinalASM.Diploid.kmer{kmer_size}"
        
    curr_time = int(time.time())
    run(
        kmer_size,
        path_fasta,
        path_save_prefix,
        f"{dir_workspace}/Flagger_Evaluation/Winnowmap_Meryl_Binning/tmp",
        f"{dir_workspace}/Flagger_Evaluation/Winnowmap_Meryl_Binning/sh",
        f"{dir_workspace}/Flagger_Evaluation/Winnowmap_Meryl_Binning/log",
        f"Meryl_Binning.KOREF_FinalASM.Diploid.kmer{kmer_size}.{curr_time}",
        96,
        memory = 3000,
        hostname = [],
        jobs_wait = []
    )