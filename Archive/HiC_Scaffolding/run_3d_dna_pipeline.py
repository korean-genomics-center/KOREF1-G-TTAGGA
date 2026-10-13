#%%
import subprocess, os, time

PATH_3DDNAPIPELINE = "/BiO/Research/Project2/Endangered_Goral/Resources/Tools/3D_DNA/3D_DNA_v201008/3d-dna-201008/run-asm-pipeline.sh"
DIR_GAWK = "/BiO/Access/yoonsung/anaconda3/bin"


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

def write_shell_script_and_run(path_shell, shell_context, env = None):
    with open(path_shell, 'w') as fw:
        fw.write(shell_context)
    
    subprocess.run(f"qsub {path_shell}", shell = True, env = env)  

def make_symlink(path_src, path_dst):
    if os.path.exists(path_dst):
        print(f"Already Exists!: {path_dst}\nSkip Symlink...")
    else:
        os.symlink(path_src, path_dst)

def run(path_asm, path_mnd, dir_save, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):
    os.makedirs(dir_save, exist_ok=True)
    
    os.chdir(dir_save)

    cmd = f"""cd {dir_save};{PATH_3DDNAPIPELINE} {path_asm} {path_mnd}"""
    
    os.makedirs(dir_tmp, exist_ok=True)
    os.makedirs(dir_shell, exist_ok=True)
    os.makedirs(dir_log, exist_ok=True)
    
    envs = os.environ
    envs["PATH"] = f'{envs["PATH"]}:{DIR_GAWK}'
    
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
        ),
        env = envs
    )

if __name__ == "__main__":
    dir_workspace = "/BiO/Research/Project2/Endangered_Goral/Analysis"
    
    for hap in [1, 2]:
        path_asm = f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Initial_Assembly/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.asm.hic.hap{hap}.p_ctg.fasta"
        
        path_mnd = f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/YaHS_Scaffolding/HiC_Alignment.Juicer/Hap{hap}/aligned/merged_nodups.txt"
        
        dir_save = f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/3D_DNA_Scaffolding/Hap{hap}"
        
        curr_time = int(time.time())
        run(
            path_asm,
            path_mnd,
            dir_save,
            f"{dir_workspace}/3D_DNA_Scaffolding/tmp",
            f"{dir_workspace}/3D_DNA_Scaffolding/sh",
            f"{dir_workspace}/3D_DNA_Scaffolding/log",
            f"Job.3D_DNA_Scaffolding.KoGoral.DraftASM.Hap{hap}.{curr_time}",
            128,
            [],
            []
        )
#%%