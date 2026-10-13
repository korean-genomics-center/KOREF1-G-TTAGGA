import os, time
import subprocess

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

def run(path_input_bam, path_output_bam, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list()):    

    cmd = f"""{PATH_SAMTOOLS} sort -@ {n_jobs} -o {path_output_bam} {path_input_bam}; {PATH_SAMTOOLS} index -@ {n_jobs} {path_output_bam}"""
        
    os.makedirs(os.path.dirname(path_output_bam), exist_ok=True)
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

        path_input_bam = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/CRAQ_Evaluation/Align.LongRead.Winnowmap/Winnowmap_Align.Ref.KOREF_FinalASM.hap{hap}.QRY.ONT_Hifi_Merged.bam"

        path_output_bam = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/CRAQ_Evaluation/Align.LongRead.Winnowmap/Winnowmap_Align.Ref.KOREF_FinalASM.hap{hap}.QRY.ONT_Hifi_Merged.sorted.bam"
        # path_output_bam = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Read_Alignment/Align.Winnowmap2/Align_PolishedASM.Polish_R8.SV/Primary_Alignment/Sorted_Bam/Winnowmap_Align.Primary.Ref.KOREF_PolishedASM.Polish_R8.SV.hap{hap}.QRY.{tgs_name}.triobin_split.triobin_{triobins}.sorted.bam"

        assert os.path.exists(path_input_bam)
        run(
            path_input_bam, 
            path_output_bam, 
            f"{dir_workspace}/CRAQ_Evaluation/Sort_n_Index_BAM/tmp", 
            f"{dir_workspace}/CRAQ_Evaluation/Sort_n_Index_BAM/sh", 
            f"{dir_workspace}/CRAQ_Evaluation/Sort_n_Index_BAM/log", 
            f"Sort_n_Index_BAM.Ref.KOREF_FinalASM.hap{hap}.ONT_Hifi_Merged.{curr_time}",
            64,
            jobs_wait=[]
        )