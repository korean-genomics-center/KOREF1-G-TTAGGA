import pandas as pd
import os, time
import subprocess

PATH_STAR = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/STAR-2.7.10b/bin/Linux_x86_64/STAR"


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

def run(dir_reference, path_read1, path_read2, path_output_prefix, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):
    cmd_filtervariant = f"""{PATH_STAR} --runThreadN {n_jobs} --twopassMode Basic --outSAMstrandField None --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {path_output_prefix} --genomeDir {dir_reference} --genomeLoad NoSharedMemory --readFilesCommand zcat --readFilesIn {path_read1} {path_read2}"""
    
    for param_name, param_value in additional_params.items():
        cmd_filtervariant += f" --{param_name} {param_value}"
    
    os.makedirs(os.path.dirname(path_output_prefix), exist_ok=True)
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
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis"
    curr_time = int(time.time())
    
    
    dict_dataprefix_format = {
        "Blood":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WTS/RemoveAdapter/Blood/KOREF-S1_Blood_WTS_NovaSeq_20190403_L01",
        "CellLine":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WTS/RemoveAdapter/CellLine/KOREF-S1_CellLine_WTS_HiSeq_20180417_L01"
    }

    for datatype, path_read_prefix in dict_dataprefix_format.items():
        for hap in [1, 2]:
            
            dir_reference = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/ShortRead_Alignment/STAR_Reference/Hap{hap}"

            path_read1 = f"{path_read_prefix}_R1.trimPolyG.Trim_F15_T5.fq.gz"
            path_read2 = f"{path_read_prefix}_R2.trimPolyG.Trim_F15_T5.fq.gz"

            path_output_prefix = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/ShortRead_Alignment/STAR_Alignment/Hap{hap}/STAR_Align.Ref_KOREF1_FinalASM_Hap{hap}.QRY_{datatype}_Trim_F15_T5.twopassMode"

            curr_time = int(time.time())
            run(
                dir_reference, path_read1, path_read2, path_output_prefix,
                f"{dir_workspace}/CAT_Annotation/STAR_Alignment/tmp",
                f"{dir_workspace}/CAT_Annotation/STAR_Alignment/sh",
                f"{dir_workspace}/CAT_Annotation/STAR_Alignment/log",
                f"STAR_Alignment.KOREF_FinalASM.Hap{hap}.{datatype}.{curr_time}",
                16,
                [],
                []
            )