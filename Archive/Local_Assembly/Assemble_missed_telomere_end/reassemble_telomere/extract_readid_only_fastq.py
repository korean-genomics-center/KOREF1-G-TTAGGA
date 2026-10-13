from pathlib import Path
import pandas as pd

import os, time
import subprocess


PATH_SCRIPT = "/BiO/Research/Project2/Endangered_Goral/Resources/Scripts/KOGIC_KOREF_AGGA/Local_Assembly/get_subseq_of_fastq.py"

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

def run(path_fastq, path_readid, path_save, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list()):    

    cmd = f"""python {PATH_SCRIPT} --input {path_fastq} --file {path_readid} --out {path_save}"""
        
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
            hostname=hostname,
            hold_jids=jobs_wait
        )
    )

if __name__ == "__main__":
    dir_workspace = "/BiO/Research/Project2/Endangered_Goral/Analysis"
      
    
    path_reads_format = "/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/Collect_ONT_with_telomere/KoGoral.Merged.ONT_reads.Adapt_Filt.Edge_Trim_Front_500_Tail_500.Filter_L50kb_QV15.telomere_count.highrepeat_500.readid.unaligned_reads.ratio_{ratio}.fastq"
    
    dict_ratio_to_comb = {
        "60":[
            "k17.min_9.max_42"
        ],
        "65":[
            "k17.min_12.max_26",
            "k17.min_12.max_42",
            "k21.min_12.max_19",
            "k21.min_12.max_26"
        ]
    }
    
    curr_time = int(time.time())
    
    for ratio in ["60", "65"]:
        for comb in dict_ratio_to_comb[ratio]:
            
            table_classified = pd.read_csv(f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/Classify_ONT_reads/KoGoral.Merged.ONT_reads.Adapt_Filt.Edge_Trim_Front_500_Tail_500.Filter_L50kb_QV15.telomere_count.highrepeat_500.readid.unaligned_reads.ratio_{ratio}.{comb}.count.tsne", sep = '\t')
            
            for cluster in table_classified["Cluster"].unique():
                if cluster == -1:
                    continue
                
                path_save_readids = f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/Clustered_ONT_Reads/KoGoral.Merged.ONT_reads.Edge_Trim_Front_500_Tail_500.Filter_L50kb_QV15.telomere_count.highrepeat_500.unaligned_reads.ratio_{ratio}.{comb}.cluster_{cluster}.readid.txt"
                
                with open(path_save_readids, 'w') as fw:
                    for readid in table_classified[table_classified["Cluster"] == cluster]["readid"]:
                        fw.write(f"{readid}\n")
            
                path_save_fastq = f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/Clustered_ONT_Reads/KoGoral.Merged.ONT_reads.Edge_Trim_Front_500_Tail_500.Filter_L50kb_QV15.telomere_count.highrepeat_500.unaligned_reads.ratio_{ratio}.{comb}.cluster_{cluster}.readid.fastq"
                
                run(
                    path_reads_format.format(ratio = ratio), 
                    path_save_readids, 
                    path_save_fastq,
                    f"{dir_workspace}/Extract_Fastq/tmp",
                    f"{dir_workspace}/Extract_Fastq/sh",
                    f"{dir_workspace}/Extract_Fastq/log",
                    f"Extract_Fastq.ONT_reads.Filter_L50kb_QV15.ratio_{ratio}.{comb}.cluster_{cluster}.{curr_time}",
                    1,
                    hostname = [],
                    jobs_wait = []
                )
