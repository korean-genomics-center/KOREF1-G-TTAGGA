import os, time, math
import subprocess
from copy import deepcopy

DIR_HMMER = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/HumASHMMER/Prerequisite/hmmer-3.4/bin/bin"
DIR_BEDOPS = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/HumASHMMER/Prerequisite/bedops_v2.4.42"
DIR_HUMAS_HMMER = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/HumASHMMER/HumASHMMER_commitf7e1e13/HumAS-HMMER_for_AnVIL"


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

def write_shell_script_and_run(path_shell, shell_context, env = None):
    with open(path_shell, 'w') as fw:
        fw.write(shell_context)
    subprocess.run(f"qsub {path_shell}", shell = True, env = env)

def make_symlink(path_src, dir_dst, ch_extension = None):
    filename_dst = os.path.basename(path_src)
    if ch_extension:
        filename_dst = f"{'.'.join(filename_dst.split('.')[:-1])}.{ch_extension}"
    path_dst = os.path.join(dir_dst, filename_dst)
    if os.path.exists(path_dst):
        pass
    else:
        os.symlink(path_src, path_dst)
    return filename_dst

def copy_data(path_src, dir_dst, ch_extension = None):
    filename_dst = os.path.basename(path_src)
    if ch_extension:
        filename_dst = f"{'.'.join(filename_dst.split('.')[:-1])}.{ch_extension}"
    path_dst = os.path.join(dir_dst, filename_dst)
    if os.path.exists(path_dst):
        pass
    else:
        subprocess.run(f"cp {path_src} {path_dst}", shell = True)
    return filename_dst

def run(path_asm, dir_out, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list()):
    os.makedirs(dir_out, exist_ok=True)
    
    copy_data(path_asm, dir_out, ch_extension="fa")
    
    list_filenames_copy = [
        "AS-HORs-hmmer3.4-071024.hmm",
        "AS-SFs-hmmer3.0.290621.hmm",
        "hmmertblout2bed.awk",
        "overlap_filter.py"
    ]
    
    for filename in list_filenames_copy:
        copy_data(os.path.join(DIR_HUMAS_HMMER, filename), dir_out)

    cmd = f"""
cd {dir_out};
{DIR_HUMAS_HMMER}/hmmer-run.sh {dir_out} AS-HORs-hmmer3.4-071024.hmm {n_jobs};
{DIR_HUMAS_HMMER}/hmmer-run_SF.sh {dir_out} AS-SFs-hmmer3.0.290621.hmm {n_jobs};
"""
    humas_env = os.environ
    humas_env["PATH"] = f"{humas_env['PATH']}:{DIR_BEDOPS}:{DIR_HMMER}"
    
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

    curr_time = int(time.time())

    
    for hap in [1, 2]:
        path_asm = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Final_Assembly/FinalASM.KOREF_S1.Trio_Specific.hap{hap}.Reordered.fasta"
        dir_out = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Repeat_Annotation/Centromere_Annotation/HumAS_HMMER/Hap{hap}"
        
        jobname = f"KOREF_ASM_Hap{hap}"
        
        run(
            path_asm, dir_out, 
            f"{dir_workspace}/Repeat_Annotation/HumAS_HMMER_Annotation/tmp",
            f"{dir_workspace}/Repeat_Annotation/HumAS_HMMER_Annotation/sh",
            f"{dir_workspace}/Repeat_Annotation/HumAS_HMMER_Annotation/log",
            f"HumAS_HMMER_Annotation.{jobname}.{curr_time}",
            64,
            [],
            []
        )
#%%