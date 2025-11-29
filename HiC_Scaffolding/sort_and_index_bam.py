import pandas as pd
import os, time
import subprocess

PATH_SAMTOOLS = "/BiO/Research/Project2/Endangered_Jindo_Dog/Resources/Tools/samtools/samtools_v1.17/samtools-1.17/bin/samtools"


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
    dir_workspace = "/BiO/Research/Project2/Endangered_Goral/Analysis"

    dict_hap_to_parental = {
        1:"p",
        2:"m"
    }
    
    dict_hap_exchange = {
        'p' : 'm',
        'm' : 'p'
    }
    
    curr_time = int(time.time())

    for hap in [1, 2]:
        path_input_bam = f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/YaHS_Scaffolding/HiC_Alignment.4Enz/Hap{hap}/results/bowtie_results/bwt2/sample/KoGoral.Merged.HiC_reads.fastp.trim_F15_T5._KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.asm.hic.hap{hap}.p_ctg.bwt2pairs.bam"

        path_output_bam = f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/YaHS_Scaffolding/HiC_Alignment.4Enz/Hap{hap}/results/bowtie_results/bwt2/sample/KoGoral.Merged.HiC_reads.fastp.trim_F15_T5._KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.asm.hic.hap{hap}.p_ctg.bwt2pairs.sorted.bam"

        assert os.path.exists(path_input_bam)
        run(
            path_input_bam, 
            path_output_bam, 
            f"{dir_workspace}/Run_YaHS/Sort_n_Index_BAM/tmp", 
            f"{dir_workspace}/Run_YaHS/Sort_n_Index_BAM/sh", 
            f"{dir_workspace}/Run_YaHS/Sort_n_Index_BAM/log", 
            f"Sort_n_Index_BAM.Hic_Alignment.hap{hap}.{curr_time}",
            48,
            jobs_wait=[],
            hostname = ["Abe", "Crane"]
        )