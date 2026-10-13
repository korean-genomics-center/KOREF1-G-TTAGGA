#%%
#%%
import pandas as pd
import os, time
import subprocess

PATH_HICCONVERTFORMAT = "/BiO/Access/yoonsung/anaconda3/envs/hicexplorer/bin/hicConvertFormat"


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

def run(path_matrice, path_output, path_hicpro_bed, dir_tmp, dir_shell, dir_log, jobname, n_jobs, outputformat = "h5", hostname = list(), jobs_wait = list()):    

    cmd = f"""{PATH_HICCONVERTFORMAT} --matrices {path_matrice} --outFileName {path_output} --inputFormat hicpro --outputFormat {outputformat} --bedFileHicpro {path_hicpro_bed}"""
        
    os.makedirs(os.path.dirname(path_output), exist_ok=True)
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

    for window in [20000, 40000, 150000, 500000, 1000000]:
        for hap in [1, 2]:
            path_matrice = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/HiC_Evaluation/HiC_Alignment/NovaSeq_20250113_Arima.MapQ0.4Enz/Hap{hap}/results/hic_results/matrix/sample/raw/{window}/sample_{window}.matrix"
            path_hicpro_bed = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/HiC_Evaluation/HiC_Alignment/NovaSeq_20250113_Arima.MapQ0.4Enz/Hap{hap}/results/hic_results/matrix/sample/raw/{window}/sample_{window}_abs.bed"

            path_output = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/HiC_Evaluation/HiC_Alignment/NovaSeq_20250113_Arima.MapQ0.4Enz/Hap{hap}/HiCPro2HiCExplorer/Matrix/sample_{window}.HiCExplorer.h5"

            assert os.path.exists(path_matrice)
            assert os.path.exists(path_hicpro_bed)

            curr_time = int(time.time())
            run(
                path_matrice, path_output, path_hicpro_bed,
                f"{dir_workspace}/HiC_Evaluation/Convert_HicPro2HiCExplorer/tmp",
                f"{dir_workspace}/HiC_Evaluation/Convert_HicPro2HiCExplorer/sh",
                f"{dir_workspace}/HiC_Evaluation/Convert_HicPro2HiCExplorer/log",
                f"Convert_HicPro2HiCExplorer.REF.KOREF_FinalASM.hap{hap}.NovaSeq_20250113_Arima.MapQ0.window_{window}.{curr_time}",
                4,
                outputformat = "h5",
                hostname = [],
                jobs_wait = []
            )