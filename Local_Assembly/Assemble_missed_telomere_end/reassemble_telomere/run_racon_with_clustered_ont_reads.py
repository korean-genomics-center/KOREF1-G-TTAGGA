#%%
import os, time
import subprocess

PATH_MINIMAP2 = "/BiO/Research/Project2/Endangered_Jindo_Dog/Resources/Tools/minimap2/minimap2_v2.30_r1287/minimap2/minimap2"
PATH_RACON = "/BiO/Research/Project2/Endangered_Jindo_Dog/Resources/Tools/Racon/Racon_v1.4.3/racon/build/bin/racon"


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

def run(other_fastq, longest_fasta, path_save_prefix, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):
    os.makedirs(os.path.dirname(path_save_prefix), exist_ok=True)

    cmd = f"""
# {PATH_MINIMAP2} -x map-ont -c -w 5 -U 1,1000000 -e 0 -t {n_jobs} {longest_fasta} {other_fastq} > {path_save_prefix}.others_align.paf; 
{PATH_RACON} --include-unpolished --no-trimming --threads {n_jobs} {other_fastq} {path_save_prefix}.others_align.paf {longest_fasta} > {path_save_prefix}.racon.fa
"""
    
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

#%%  
if __name__ == "__main__":
    dir_workspace = "/BiO/Research/Project2/Endangered_Goral/Analysis"

    curr_time = int(time.time())

    for ind in range(7):
        
        longest_fasta = f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/Cluster_Unaligned_ONT_Reads_by_Alignment/ONT_Reads.Telomeric_Read_Self_bothhap_a2h80idy95_ratio_60.mapontPcw5U11000000e0.idy95align1000hang100000a2h90.lenover_100000.alignment_ratio_90.readset_{ind}.longest.fa"
        
        other_fastq = f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/Cluster_Unaligned_ONT_Reads_by_Alignment/ONT_Reads.Telomeric_Read_Self_bothhap_a2h80idy95_ratio_60.mapontPcw5U11000000e0.idy95align1000hang100000a2h90.lenover_100000.alignment_ratio_90.readset_{ind}.others.fastq"
        
        jobname = f"ont_reads.cluster_{ind}"

        path_save_prefix = f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/Consensus_read/Racon_Consensus.Telomeric_Read_Self_bothhap_a2h80idy95_ratio_60.idy95align1000hang100000a2h90.lenover_100000.alignment_ratio_90.readset_{ind}"

        assert os.path.exists(longest_fasta)
        assert os.path.exists(other_fastq)
        
        run(
            other_fastq, longest_fasta, path_save_prefix, 
            f"{dir_workspace}/RACON_Consensus/tmp",
            f"{dir_workspace}/RACON_Consensus/sh",
            f"{dir_workspace}/RACON_Consensus/log",
            f"RACON_Consensus.{jobname}.{curr_time}",
            128,
            [],
            []
        )