#%%
import subprocess, os, time

PATH_YAHS = "/BiO/Research/Project2/Endangered_Goral/Resources/Tools/YaHS/YaHS-v1.2.2/yahs/yahs"


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

def run(path_contig, path_bam, path_output_prefix, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):

    cmd = f"{PATH_YAHS} -o {path_output_prefix} --no-contig-ec --no-scaffold-ec --file-type BAM {path_contig} {path_bam}"
    
    os.makedirs(os.path.dirname(path_output_prefix), exist_ok=True)
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
    dir_workspace = "/BiO/Research/Project2/Endangered_Goral/Analysis"
    
    for hap in [1]:
        path_contig = f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/HiC_Scaffolding_with_additional_contigs/HiC_Alignment.4Enz.with_raw_unclustered_read/genome/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.asm.hic.hap1.p_ctg.with_additional_teloseq.with_unclustered_rawONT.fasta"
        path_bam = f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/HiC_Scaffolding_with_additional_contigs/HiC_Alignment.4Enz.with_raw_unclustered_read/results/bowtie_results/bwt2/sample/KoGoral.Merged.HiC_reads.fastp.trim_F15_T5._KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.asm.hic.hap1.p_ctg.with_additional_teloseq.with_unclustered_rawONT.bwt2pairs.Sorted.MarkDuplicates.SortByName.bam"
        
        path_output_prefix = f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/HiC_Scaffolding_with_additional_contigs/YaHS_Run.with_additional_teloseq.with_raw_unclustered_read/YaHS_Run.KoGoral_DraftASM.Hap1.with_additional_teloseq.with_raw_unclustered_read"
        
        curr_time = int(time.time())
        run(
            path_contig,
            path_bam,
            path_output_prefix,
            f"{dir_workspace}/Run_YaHS/tmp",
            f"{dir_workspace}/Run_YaHS/sh",
            f"{dir_workspace}/Run_YaHS/log",
            f"Run_YaHS.KoGoral.DraftASM.Hap{hap}.with_additional_teloseq.with_raw_unclustered_read.{curr_time}",
            192,
            [],
            []
        )
#%%