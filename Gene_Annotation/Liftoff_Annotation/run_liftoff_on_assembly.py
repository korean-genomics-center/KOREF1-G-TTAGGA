import os, time
import subprocess
from copy import deepcopy

PATH_LIFTOFF = "/BiO/Access/yoonsung/anaconda3/envs/liftoff/bin/liftoff"
PATH_MINIMAP2 = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/minimap2/minimap2_2.29_r1283/minimap2-2.29_x64-linux/minimap2"


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

def make_symlink(path_src, dir_dst, ch_extension = None):
    filename_dst = os.path.basename(path_src)
    if ch_extension:
        filename_dst = f"{'.'.join(filename_dst.split('.')[:-1])}.{ch_extension}"
    path_dst = os.path.join(dir_dst, filename_dst)
    if os.path.exists(path_dst):
        pass
    else:
        os.symlink(path_src, path_dst)
    return filename_dst

def run(path_asm_qry, path_ref, path_gff_db, dict_chromosomes_qry_to_ref, path_output_prefix, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):
    dir_save = os.path.dirname(path_output_prefix)
    os.makedirs(dir_save, exist_ok=True)
    
    path_chroms = f"{path_output_prefix}.chroms.txt"
    path_save_output_gff = f"{path_output_prefix}.annotate.gff3"
    path_save_unmapped = f"{path_output_prefix}.unmapped"
    dir_save_intermediate = os.path.join(os.path.dirname(path_output_prefix), "intermediate")
    
    with open(path_chroms, 'w') as fw:
        for chr_qry, chr_ref in dict_chromosomes_qry_to_ref.items():
            fw.write(f"{chr_qry},{chr_ref}\n")

    cmd = f"{PATH_LIFTOFF} {path_asm_qry} {path_ref} -db {path_gff_db} -o {path_save_output_gff} -u {path_save_unmapped} -dir {dir_save_intermediate} -chroms {path_chroms} -p {n_jobs} -m {PATH_MINIMAP2} -sc 0.95 -copies -polish"
    
    for key, val in additional_params.items():
        cmd += f"--{key} {val} "
    
    os.makedirs(dir_save_intermediate, exist_ok=True)
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

    modified_option = {
    }
    
    path_ref = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/CHM13/chm13v2.0.fa"
    path_gff_db = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/CHM13/Gene_Annotation/chm13.draft_v2.0.gene_annotation.gff3_db"
    
    list_chr_order = list(map(lambda val: f"chr{val}", range(1, 23))) # + ["chrX", "chrY"]
    
    for hap in [1, 2]:
        path_asm_qry = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Final_Assembly/FinalASM.KOREF_S1.Trio_Specific.hap{hap}.Reordered.fasta"
        
        path_output_prefix = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/Liftoff/From.chm13v2_CATLiftoff_v2.To.KOREF_ASM_20251027_Hap{hap}/liftoff_annotation.KOREF_ASM_Hap{hap}.from_chm13v2_CATLiftoff_v2"
        
        dict_chromosomes_qry_to_ref = dict(zip(list_chr_order, list_chr_order))
        if hap == 1:
            dict_chromosomes_qry_to_ref["chrY"] = "chrY"
        else:
            dict_chromosomes_qry_to_ref["chrX"] = "chrX"
    
        jobname = f"Liftoff_KOREF_ASM_Hap{hap}.from_chm13v2_CATLiftoff_v2"
        
        run(
            path_asm_qry,
            path_ref,
            path_gff_db,
            dict_chromosomes_qry_to_ref,
            path_output_prefix, 
            f"{dir_workspace}/Gene_Annotation/Liftoff_Annotation/tmp",
            f"{dir_workspace}/Gene_Annotation/Liftoff_Annotation/sh",
            f"{dir_workspace}/Gene_Annotation/Liftoff_Annotation/log",
            f"Liftoff_Annotation.{jobname}.{curr_time}",
            12,
            [],
            [],
            **modified_option
        )
#%%