#%%
import subprocess, os, time

PATH_CACTUSPAN = "/BiO/Access/yoonsung/anaconda3/envs/cactus/bin/cactus-pangenome"

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
    
def run_job(path_reference, path_assembly, name_ref, name_asm, dir_save, name_out, jobsuffix, dir_tmp, dir_shell, dir_log, n_jobs = 1):
    dir_jobscript = os.path.join(dir_save, f"js.{int(time.time())}")
    
    os.makedirs(dir_save, exist_ok=True)
    
    path_save_seqfile = os.path.join(dir_save, f"{name_out}.seqFile.txt")
    with open(path_save_seqfile, 'w') as fw:
        fw.write(f"{name_ref}\t{path_reference}\n")
        fw.write(f"{name_asm}\t{path_assembly}\n")
    
    cmd = f"""
{PATH_CACTUSPAN} {dir_jobscript} {path_save_seqfile} --outDir {dir_save} --outName {name_out} --reference {name_ref} --maxCores {n_jobs} --mgCores {n_jobs} --mapCores {max(1, int(n_jobs/4))} --consCores {max(1, int(n_jobs/4))} --mgCores {n_jobs} --indexCores {max(1, int(n_jobs/2))} 
"""
    
    jobname = f"Run_Cactus_Pangenome_{jobsuffix}.{int(time.time())}"
    
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
    workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis/CAT_Annotation/Cactus_Alignment"
    
    path_reference = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Cactus_Alignment/GCF_000001405.40_GRCh38.p14_genomic.major_chromosome.id_conv.fa"
    name_ref = "GRCh38"
    
    for hap in [1, 2]:
        path_assembly = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Repeat_Annotation/RepeatMasker/Hap{hap}/FinalASM.KOREF_S1.Trio_Specific.hap{hap}.Reordered.fasta.masked"
        name_asm = f"KOREF1_ASM_Hap{hap}"
        
        dir_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Cactus_Alignment/Align.Hap{hap}"
        name_out = f"grch38_koref1hap{hap}"
        
        run_job(
            path_reference, path_assembly, name_ref, name_asm, dir_save, name_out,
            f"KOREF1_G_ASM.Hap{hap}.Masked",
            os.path.join(workspace, "tmp"),
            os.path.join(workspace, "sh"),
            os.path.join(workspace, "log"),
            n_jobs = 128
        )