import pandas as pd
import os, time
import subprocess

PATH_FLYE = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/Flye/Flye_v2.9.6/Flye/bin/flye"


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

def run(path_hifi, dir_out, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **modified_arguments):    

    cmd = f"""{PATH_FLYE} --pacbio-hifi {path_hifi} -o {dir_out} -t {n_jobs}"""
    
    for key, val in modified_arguments.items():
        if len(key) == 1:
            cmd += f" -{key} {val}"
        else:
            cmd += f" --{key} {val}"
        
    os.makedirs(dir_out, exist_ok=True)
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
    
    dict_hap_to_parental = {
        1:'p',
        2:'m'
    }

    for hap in dict_hap_to_parental.keys():
        dir_out = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20.with_0/Hap{hap}"
            
        path_hifi = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/Input_Manipulation/Parental_Split_Data/PacBio_HiFi_ccs_converted.filt.Filter_L1000_QV20.merged.triobin_split.triobin_{dict_hap_to_parental[hap]}.fq.gz"
        path_hifi_0 = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/Input_Manipulation/Parental_Split_Data/PacBio_HiFi_ccs_converted.filt.Filter_L1000_QV20.merged.triobin_split.triobin_0.fq.gz"

        assert os.path.exists(path_hifi)
        assert os.path.exists(path_hifi_0)

        curr_time = int(time.time())
        run(
            f"{path_hifi} {path_hifi_0}", 
            dir_out,
            f"{dir_workspace}/Flye_Binned_Assembly/tmp",
            f"{dir_workspace}/Flye_Binned_Assembly/sh",
            f"{dir_workspace}/Flye_Binned_Assembly/log",
            f"Flye_Binned_Assembly.HiFi.Filter_L1000_QV20.hap{hap}.with_0.{curr_time}",
            64,
            hostname = [],
            jobs_wait = []
        )
