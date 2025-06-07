#%%
import subprocess, os, time

PATH_SEQTK = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/SeqTK/seqtk/seqtk"

def get_context_for_qsub_submission(cmd, jobname, stderr_path, stdout_path, nthreads, hold_jids = list()):
    #Header
    script_header = f"#!bin/sh\n"
    script_header += f"#$ -N {jobname}\n"
    script_header += f"#$ -S /bin/bash\n"
    script_header += f"#$ -V\n"
    if len(hold_jids) > 0:
        script_header += f"#$ -hold_jid {','.join(hold_jids)}\n"
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
    
def run_job(path_fastq, path_save, jobname, dir_tmp, dir_shell, dir_log, n_jobs = 1):
    cmd = f"""{PATH_SEQTK} seq -a {path_fastq} > {path_save}"""
    
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
            list()
        )
    )


if __name__ == "__main__":
    workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis"
    
    for parental in ['a', '0']:
        path_fastq = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/ONT/Parental_Split/ONT_UL.Adapt_Filt.Edge_Trim_Front_1500_Tail_1500.Filter_L50kb_QV20_Q30Perc80.merged.triobin_split.triobin_{parental}.fq.gz"
        path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/ONT/Parental_Split/ONT_UL.Adapt_Filt.Edge_Trim_Front_1500_Tail_1500.Filter_L50kb_QV20_Q30Perc80.merged.triobin_split.triobin_{parental}.fasta"
        
        run_job(
            path_fastq,
            path_save,
            f"Convert_FQ2FA.ONT_UL.Edge_Trim_Front_1500_Tail_1500.Filter_L50kb_QV20_Q30Perc80.merged.triobin_split.triobin_{parental}",
            os.path.join(workspace, "ConvertFQ2FA", "tmp"),
            os.path.join(workspace, "ConvertFQ2FA", "sh"),
            os.path.join(workspace, "ConvertFQ2FA", "log"),
            1
        )
    
    for parental in ['a', '0']:
        path_fastq = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/PacBio/Parental_Split/PacBio_HiFi_ccs_converted.filt.Filter_L1000_QV20.merged.triobin_split.triobin_{parental}.fq.gz"
        path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/PacBio/Parental_Split/PacBio_HiFi_ccs_converted.filt.Filter_L1000_QV20.merged.triobin_split.triobin_{parental}.fasta"
        
        run_job(
            path_fastq,
            path_save,
            f"Convert_FQ2FA.PacBio_HiFi_ccs_converted.Filter_L1000_QV20.merged.triobin_split.triobin_{parental}",
            os.path.join(workspace, "ConvertFQ2FA", "tmp"),
            os.path.join(workspace, "ConvertFQ2FA", "sh"),
            os.path.join(workspace, "ConvertFQ2FA", "log"),
            1
        )
#%%