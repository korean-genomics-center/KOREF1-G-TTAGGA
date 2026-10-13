import os, time
import subprocess
from copy import deepcopy

PATH_REPEATMODELER = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/RepeatMasker/RepeatModeler_v2.0.7/RepeatModeler-2.0.7/RepeatModeler"


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

def run(path_db, path_save_log, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):

    cmd = f"cd {os.path.dirname(path_db)}; {PATH_REPEATMODELER} -database {path_db} -threads {n_jobs} -LTRStruct -srand 0"
    
    for key, val in additional_params.items():
        cmd += f" {key} {val} "
    
    cmd += f" > {path_save_log} 2>&1"
    
    dir_save = os.path.dirname(path_save_log)
    os.makedirs(dir_save, exist_ok=True)
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

    modified_option = {
    }
    
    dict_hap_to_trio = {
        1:"Paternal",
        2:"Maternal"
    }
    
    for hap in [1, 2]:
        path_db = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Repeat_Annotation/RepeatModeler/Hap{hap}/DB/KOREF1_{dict_hap_to_trio[hap]}_v1"
        path_save_log = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Repeat_Annotation/RepeatModeler/Hap{hap}/DB/KOREF1_{dict_hap_to_trio[hap]}_v1.run_repeatmodeler.rerun_repeatmasker.log"
    
        jobname = f"RepeatModeler_KOREF_ASM_Hap{hap}"
        
        run(
            path_db,
            path_save_log, 
            f"{dir_workspace}/Repeat_Annotation/RepeatModeler_Annotation/tmp",
            f"{dir_workspace}/Repeat_Annotation/RepeatModeler_Annotation/sh",
            f"{dir_workspace}/Repeat_Annotation/RepeatModeler_Annotation/log",
            f"RepeatModeler_Annotation.{jobname}.{curr_time}",
            48,
            [],
            [],
            **modified_option
        )
#%%