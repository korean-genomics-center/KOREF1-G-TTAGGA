#%%
import subprocess, os, time

PATH_LIMA = "/BiO/Access/yoonsung/anaconda3/envs/pacbio/bin/lima"

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
    
def run_job(path_bam, path_primer, path_save, jobsuffix, dir_tmp, dir_shell, dir_log, n_jobs = 1):
    
    cmd = f"""
{PATH_LIMA} --isoseq --peek-guess -j {n_jobs} {path_bam} {path_primer} {path_save}
"""
    
    jobname = f"Run_Lima_{jobsuffix}.{int(time.time())}"
    
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
    workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis/CAT_Annotation/Lima_on_CCS"
    
    path_primerset = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/IsoSeq_Preprocessing/Primers/Old_Isoseq_Primer.fasta"
    
    
    path_bam = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/IsoSeq_Preprocessing/Subread_to_CCS/KOREF-S1_Blood_WTS_Sequel_20181122_FC0{fc}.ccs.bam"
    path_save = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/IsoSeq_Preprocessing/Lima_isoseq/KOREF-S1_Blood_WTS_Sequel_20181122_FC0{fc}.ccs.fl.bam"
    
    
    for fc in range(5):
        run_job(
            path_bam.format(fc = fc),
            path_primerset,
            path_save.format(fc = fc),
            f"KOREF-S1_Blood_WTS_Sequel_20181122_FC0{fc}",
            os.path.join(workspace, "tmp"),
            os.path.join(workspace, "sh"),
            os.path.join(workspace, "log"),
            n_jobs = 8
        )
    
    path_bam = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/IsoSeq_Preprocessing/Subread_to_CCS/KOREF-S1_CellLine_WTS_Sequel_20181122_FC0{fc}.ccs.bam"
    path_save = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/IsoSeq_Preprocessing/Lima_isoseq/KOREF-S1_CellLine_WTS_Sequel_20181122_FC0{fc}.ccs.fl.bam"
    
    
    for fc in range(5):
        run_job(
            path_bam.format(fc = fc),
            path_primerset,
            path_save.format(fc = fc),
            f"KOREF-S1_CellLine_WTS_Sequel_20181122_FC0{fc}",
            os.path.join(workspace, "tmp"),
            os.path.join(workspace, "sh"),
            os.path.join(workspace, "log"),
            n_jobs = 8
        )
    
    path_bam = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/IsoSeq_Preprocessing/Subread_to_CCS/KOREF-S1_Blood_WTS_Sequel-II_20230519_FC00.ccs.bam"
    path_save = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/IsoSeq_Preprocessing/Lima_isoseq/KOREF-S1_Blood_WTS_Sequel-II_20230519_FC00.ccs.fl.bam"
    run_job(
        path_bam,
        path_primerset,
        path_save,
        f"KOREF-S1_Blood_WTS_Sequel-II_20230519_FC00",
        os.path.join(workspace, "tmp"),
        os.path.join(workspace, "sh"),
        os.path.join(workspace, "log"),
        n_jobs = 8
    )