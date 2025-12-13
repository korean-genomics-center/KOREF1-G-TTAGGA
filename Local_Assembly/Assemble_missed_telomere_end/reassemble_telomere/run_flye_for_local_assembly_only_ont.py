import pandas as pd
import os, time
import subprocess

from glob import glob

PATH_FLYE = "/BiO/Research/Project2/Endangered_Jindo_Dog/Resources/Tools/flye/flye_v2.9.6/Flye/bin/flye"


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

def run(path_ont_hq, dir_out, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):    

    cmd = f"""{PATH_FLYE} --nano-hq {path_ont_hq} -o {dir_out} -t {n_jobs}"""
    
    for key, val in additional_params.items():
        cmd += f" {key} {val}"
    
    os.makedirs(dir_out, exist_ok=True)
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
    
    
    for align_ratio in [75, 80]:
        path_ont_src = f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/Collect_ONT_with_telomere/KoGoral.Merged.ONT_reads.Adapt_Filt.Edge_Trim_Front_500_Tail_500.Filter_L50kb_QV15.telomere_count.highrepeat_500.readid.unaligned_reads.ratio_{align_ratio}.fastq"

        dir_out = f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/Flye_Assembly/Flye.ONT.telorepeat_500.unaligned_ratio_{align_ratio}.g300k"

        
        assert os.path.exists(path_ont_src)
        # continue

        curr_time = int(time.time())
        run(
            path_ont_src,  
            dir_out,
            f"{dir_workspace}/Flye_Local_Assembly/tmp",
            f"{dir_workspace}/Flye_Local_Assembly/sh",
            f"{dir_workspace}/Flye_Local_Assembly/log",
            f"Flye_Local_Assembly.ONT.unaligned_{align_ratio}.{curr_time}",
            64,
            hostname = [],
            jobs_wait = [],
            **{
                "--genome-size" : "300k"
            }
        )