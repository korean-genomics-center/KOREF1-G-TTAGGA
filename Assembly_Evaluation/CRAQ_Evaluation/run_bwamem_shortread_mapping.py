#%%
import pandas as pd
import os, time
import subprocess

PATH_BWAMEM2 = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/bwa/bwa-mem_v2.3/bwa-mem2-2.3_x64-linux/bwa-mem2"
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

def run(path_fasta_ref, path_read1, path_read2, path_save, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list()):    

    cmd = f"""{PATH_BWAMEM2} mem -t {n_jobs} {path_fasta_ref} {path_read1} {path_read2} | {PATH_SAMTOOLS} view -bS -@ {min(n_jobs, 16)} -o {path_save}"""
        
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
    dict_hap_to_parental = {
        1:'p',
        2:'m'
    }
    
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis"
        
    dir_save = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/CRAQ_Evaluation/Align.ShortRead.Bwamem2"
    
    path_read1 = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/IlluminaWGS/Merged_Read/KOREF1_WGS_NovaSeq_ReadMerged_R1.trimPolyG.Trim_F15_T5.fq.gz"
    path_read2 = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/IlluminaWGS/Merged_Read/KOREF1_WGS_NovaSeq_ReadMerged_R2.trimPolyG.Trim_F15_T5.fq.gz"

    for hap in [1, 2]:
            
        path_fasta_ref = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/CRAQ_Evaluation/Align.ShortRead.Bwamem2/Reference/FinalASM.KOREF_S1.Trio_Specific.hap{hap}.Reordered.fasta"

        path_save = os.path.join(dir_save, f"Bwamem2_Align.Ref.KOREF_FinalASM.hap{hap}.QRY.IlluminaWGS.bam")

        assert os.path.exists(path_fasta_ref)
        assert os.path.exists(path_read1)
        assert os.path.exists(path_read2)

        curr_time = int(time.time())
        run(
            path_fasta_ref, path_read1, path_read2, path_save,
            f"{dir_workspace}/CRAQ_Evaluation/Bwamem2_Align/tmp",
            f"{dir_workspace}/CRAQ_Evaluation/Bwamem2_Align/sh",
            f"{dir_workspace}/CRAQ_Evaluation/Bwamem2_Align/log",
            f"Bwamem2_Align.REF.KOREF_FinalASM.hap{hap}.IlluminaWGS.{curr_time}",
            64,
            hostname = [],
            jobs_wait = []
        )