#%%
import subprocess, os, time

PATH_MINIMAP2 = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/minimap2/minimap2_2.29_r1283/minimap2-2.29_x64-linux/minimap2"
PATH_SAMTOOLS = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/Samtools/samtools-1.17/samtools-1.17/bin/samtools"

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
    
def run_job(path_reference, path_isoseq, path_save, jobsuffix, dir_tmp, dir_shell, dir_log, n_jobs = 1):
    
    cmd = f"""
{PATH_MINIMAP2} -t {n_jobs} -ax splice:hq --sam-hit-only --MD --eqx --secondary=no {path_reference} {path_isoseq} | {PATH_SAMTOOLS} view -bS -o {path_save}
"""
    
    jobname = f"Run_Minimap2_SpliceHQ_{jobsuffix}.{int(time.time())}"
    
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
    workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis/CAT_Annotation/Isoseq_Alignment"
    
    
    for hap in [1, 2]:
        for isoseqtype in ["all", "blood"]:
            path_reference = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Final_Assembly/FinalASM.KOREF_S1.Trio_Specific.hap{hap}.Reordered.fasta"
            
            path_isoseq = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/IsoSeq_Preprocessing/Isoseq_cluster/KOREF-S1.{isoseqtype}_IsoSeq.ccs.fl.refine_nc.clustered.fastq"
            
            path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/IsoSeq_Preprocessing/Isoseq_align/Minimap2_align.Ref_KOREF1_Hap{hap}.Qry_Isoseq_{isoseqtype}.splicehqMDeqxNosec.bam"
            
            run_job(
                path_reference, path_isoseq, path_save,
                f"KOREF1_G_ASM.Ref_Hap{hap}.Qry_Isoseq_{isoseqtype}",
                os.path.join(workspace, "tmp"),
                os.path.join(workspace, "sh"),
                os.path.join(workspace, "log"),
                n_jobs = 64
            )