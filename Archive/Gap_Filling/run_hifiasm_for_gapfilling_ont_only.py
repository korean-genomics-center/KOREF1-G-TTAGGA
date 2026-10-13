import pandas as pd
import os, time
import subprocess

PATH_HIFIASM = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/hifiasm_v0.25.0-r726/hifiasm/hifiasm"


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

def run(path_ont_read, path_save, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):    

    cmd = f"""{PATH_HIFIASM} --ont -o {path_save} --n-hap 1 -l0 --telo-m CCCTAA -t{n_jobs}"""
    
    for param_name, param_value in additional_params.items():
        if len(param_name) == 1:
            cmd += f" -{param_name} {param_value}"
        else:
            cmd += f" --{param_name} {param_value}"
    cmd += f" {path_ont_read}"
        
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
    
    dict_hap_to_parental = {
        1:'p',
        2:'m'
    }
    
    for hap in [2, 1]:            
        path_ont = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/Input_Manipulation/Parental_Split_Data/ONT_UL.Adapt_Filt.Edge_Trim_Front_1500_Tail_1500.Filter_L50kb_QV20_Q30Perc80.merged.triobin_split.triobin_{dict_hap_to_parental[hap]}.fq.gz"

        path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.Binned_Input/ONT_UL.Filter_L50kb/Hap{hap}/KOREF_S1.Binned_Assembly.ONT_UL.Filter_L50kb.triobin_{dict_hap_to_parental[hap]}.haploid.l0.asm"

        assert os.path.exists(path_ont), path_ont

        curr_time = int(time.time())
        run(
            path_ont,  
            path_save,
            f"{dir_workspace}/Hifiasm_Binned_Assembly/tmp",
            f"{dir_workspace}/Hifiasm_Binned_Assembly/sh",
            f"{dir_workspace}/Hifiasm_Binned_Assembly/log",
            f"Hifiasm_Binned_Assembly.ONT.hap{hap}.l0.{curr_time}",
            64,
            hostname = [],
            jobs_wait = []
        )