import pandas as pd
import os, time
import subprocess

PATH_BEDTOOLS = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/Bedtools/bedtools"


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

def run(path_censat_bed, path_output_prefix, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list()):    

    cmd = f"""
censat_patterns=("hsat1A hsat1B hsat2 hsat3 active_hor bsat")
for pattern in ${{censat_patterns[@]}};
do
    cat {path_censat_bed} | \\
        grep -i $pattern | \\
        awk '{{print $1"\\t"$2"\\t"$3}}' | \\
        {PATH_BEDTOOLS} sort -i - | \\
        {PATH_BEDTOOLS} merge -i - > {path_output_prefix}.decomposed.$pattern.bed
done

cat {path_censat_bed} | \\
    awk '$4 != "ct"{{print $1"\\t"$2"\\t"$3}}' | \\
    {PATH_BEDTOOLS} sort -i - | \\
    {PATH_BEDTOOLS} merge -i - > {path_output_prefix}.no_ct.bed
"""
        
    os.makedirs(os.path.dirname(path_output_prefix), exist_ok=True)
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

    path_censat_bed = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Flagger_Evaluation/Annotations/Suffix_Attatched_BED/CenSat.both.suffixed.bed"
    path_output_prefix = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Flagger_Evaluation/Annotations/Suffix_Attatched_BED/CenSat.both.suffixed"
    
    curr_time = int(time.time())

    run(
        path_censat_bed, 
        path_output_prefix, 
        f"{dir_workspace}/Flagger_Evaluation/Decompose_CenSat/tmp", 
        f"{dir_workspace}/Flagger_Evaluation/Decompose_CenSat/sh", 
        f"{dir_workspace}/Flagger_Evaluation/Decompose_CenSat/log", 
        f"Decompose_CenSat.KOREF1_FinalASM.{curr_time}",
        1,
        jobs_wait=[]
    )