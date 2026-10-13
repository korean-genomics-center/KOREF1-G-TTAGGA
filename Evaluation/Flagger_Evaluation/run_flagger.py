import pandas as pd
import os, time
import subprocess

PATH_WINNNOWMAP = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/winnowmap/winnowmap.v2.03_fork/Winnowmap/bin/winnowmap"
PATH_SAMTOOLS = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/Samtools/samtools-1.17/samtools-1.17/bin/samtools"


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

def run(path_coverage, path_alpha, dir_save, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list()):    

    cmd = f"""
docker run --ipc=host \\
    -v {os.path.dirname(path_coverage)}:/input \\
    -v {dir_save}:/output \\
    -v {os.path.dirname(path_alpha)}:/alphatsv \\
    mobinasri/flagger:v1.1.0 hmm_flagger \\
    --input /input/{os.path.basename(path_coverage)} \\
    --outputDir /output \\
    --trackName hmm_flagger_v1.1.0 \\
    --labelNames Err,Dup,Hap,Col \\
    --threads {n_jobs} \\
    --alphaTsv /alphatsv/{os.path.basename(path_alpha)}
"""
        
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

    dict_datatype_to_alphatsv = {
        "ONT_10kb_QV15" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Flagger_Evaluation/AlphaTsv/alpha_optimum_trunc_exp_gaussian_w_8000_n_50.ONT_R1041_Dorado_DEC_2024.v1.1.0.tsv",
        "HiFi_1kb_QV20" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Flagger_Evaluation/AlphaTsv/alpha_optimum_trunc_exp_gaussian_w_16000_n_50.HiFi_DC_1.2_DEC_2024.v1.1.0.tsv"
    }

    curr_time = int(time.time())

    for datatype in dict_datatype_to_alphatsv.keys():
        path_coverage = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Flagger_Evaluation/Bam2Cov/{datatype}/KOREF_FinalASM.Diploid.{datatype}.cov.gz"
        path_alpha = dict_datatype_to_alphatsv[datatype]
        dir_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Flagger_Evaluation/Flagger/{datatype}"
        
        run(
            path_coverage, path_alpha, dir_save,
            f"{dir_workspace}/Flagger_Evaluation/Run_Flagger/tmp",
            f"{dir_workspace}/Flagger_Evaluation/Run_Flagger/sh",
            f"{dir_workspace}/Flagger_Evaluation/Run_Flagger/log",
            f"Run_Flagger.REF.KOREF_FinalASM.Diploid.Qry.{datatype}.{curr_time}",
            64,
            hostname = [],
            jobs_wait = []
        )