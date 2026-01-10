#%%
import subprocess, os, time

PATH_ISOSEQ = "/BiO/Access/yoonsung/anaconda3/envs/pacbio/bin/isoseq"

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
    
def run_job(path_fofn, path_save, path_log, jobsuffix, dir_tmp, dir_shell, dir_log, n_jobs = 1):
    
    cmd = f"""
{PATH_ISOSEQ} cluster2 -j {n_jobs} --log-file {path_log} --log-level INFO {path_fofn} {path_save}
"""
    
    jobname = f"Run_Isoseq_cluster2_{jobsuffix}.{int(time.time())}"
    
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
    workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis/CAT_Annotation/Isoseq_cluster2"
    
    path_fofn = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/IsoSeq_Preprocessing/Isoseq_cluster/KOREF-S1.blood_IsoSeq.ccs.fl.refine_nc.sequelI.fofn"
    path_save = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/IsoSeq_Preprocessing/Isoseq_cluster/KOREF-S1.blood_IsoSeq.ccs.fl.refine_nc.sequelI.clustered.bam"
    path_log = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/IsoSeq_Preprocessing/Isoseq_cluster/KOREF-S1.blood_IsoSeq.ccs.fl.refine_nc.sequelI.clustered.log"
    
    assert os.path.exists(path_fofn)
    run_job(
        path_fofn,
        path_save,
        path_log,
        f"KOREF-S1.blood_sequelI_IsoSeq",
        os.path.join(workspace, "tmp"),
        os.path.join(workspace, "sh"),
        os.path.join(workspace, "log"),
        n_jobs = 64
    )