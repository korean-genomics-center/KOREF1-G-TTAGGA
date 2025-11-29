#%%
import subprocess, os, time

PATH_SAMTOOLS = "/BiO/Research/Project2/Endangered_Jindo_Dog/Resources/Tools/samtools/samtools_v1.17/samtools-1.17/bin/samtools"
DIR_BWA = "/BiO/Research/Project2/Endangered_Goral/Resources/Tools/bwa/bwa_0.7.19/bwa-0.7.19"
DIR_GAWK = "/BiO/Access/yoonsung/anaconda3/bin"

DIR_JUICER = "/BiO/Research/Project2/Endangered_Goral/Resources/Tools/Juicer/Juicer_v1.6_master/juicer_1.6_CPU"
PATH_GENERATESITEPOSITIONS = "/BiO/Research/Project2/Endangered_Goral/Resources/Tools/Juicer/Juicer_v1.6_master/juicer-1.6/misc/generate_site_positions.py"


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

def run(path_asm, path_hic_1, path_hic_2, dir_save, enzyme, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):
    
    dir_link_fastq = os.path.join(dir_save, "fastq")
    os.makedirs(dir_link_fastq, exist_ok=True)
    
    if path_hic_1.endswith(".gz"):
        make_symlink(path_hic_1, os.path.join(dir_link_fastq, "hic_read_R1.fastq.gz"))
        make_symlink(path_hic_2, os.path.join(dir_link_fastq, "hic_read_R2.fastq.gz"))
    else:
        make_symlink(path_hic_1, os.path.join(dir_link_fastq, "hic_read_R1.fastq"))
        make_symlink(path_hic_2, os.path.join(dir_link_fastq, "hic_read_R2.fastq"))
    
    path_copy_asm = os.path.join(dir_save, os.path.basename(path_asm))
    path_save_asm_size = os.path.join(dir_save, f"{os.path.basename(path_asm)}.chrom.sizes")
    
    cmd_process_genome = f"""
cp {path_asm} {path_copy_asm};
{PATH_SAMTOOLS} faidx {path_copy_asm};
cut -f1,2 {path_copy_asm}.fai > {path_save_asm_size};
{DIR_BWA}/bwa index {path_copy_asm};
"""

    cmd_process_ligationsite = f"""
python {PATH_GENERATESITEPOSITIONS} {enzyme} {os.path.basename(path_asm).replace(".fasta", '')} {path_copy_asm}
"""

    cmd_process_hic = f"""
{DIR_JUICER}/scripts/juicer.sh -d {dir_save} -p {path_save_asm_size} -y {path_copy_asm}_{enzyme}.txt -z {path_copy_asm} -t {n_jobs} -D {DIR_JUICER}
"""

    cmd = ""
    if not os.path.exists(f"{path_copy_asm}.bwt"):
        cmd += cmd_process_genome
    if not os.path.exists(f"{path_copy_asm.replace('.fasta', '')}_{enzyme}.txt"):
        cmd += cmd_process_ligationsite
    cmd += cmd_process_hic
    
    os.makedirs(dir_save, exist_ok=True)
    os.makedirs(dir_tmp, exist_ok=True)
    os.makedirs(dir_shell, exist_ok=True)
    os.makedirs(dir_log, exist_ok=True)
    
    envs = os.environ
    envs["PATH"] = f'{envs["PATH"]}:{DIR_BWA}:{DIR_GAWK}'
    
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
        
        path_hic_1 = "/BiO/Research/Project2/Endangered_Goral/Results/Input_Data_Filtering/ShortRead/HiC/KoGoral/Merged_Fastq/KoGoral.Merged.HiC_reads.fastp.trim_F15_T5.R1.fq.gz"
        path_hic_2 = "/BiO/Research/Project2/Endangered_Goral/Results/Input_Data_Filtering/ShortRead/HiC/KoGoral/Merged_Fastq/KoGoral.Merged.HiC_reads.fastp.trim_F15_T5.R2.fq.gz"
        
        dir_save = f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/YaHS_Scaffolding/HiC_Alignment.Juicer/Hap{hap}"
        
        enzyme="Arima4Enz"
        
        curr_time = int(time.time())
        run(
            path_asm,
            path_hic_1,
            path_hic_2,
            dir_save,
            enzyme,
            f"{dir_workspace}/HiC_Alignment/tmp",
            f"{dir_workspace}/HiC_Alignment/sh",
            f"{dir_workspace}/HiC_Alignment/log",
            f"HiC_Alignment.Juicer.KoGoral.DraftASM.Hap{hap}.{curr_time}",
            128,
            [],
            []
        )
#%%