import pandas as pd
import os, time
import subprocess

PATH_HIFIASM = "/BiO/Research/Project2/Endangered_Jindo_Dog/Resources/Tools/hifiasm/hifiasm_v0.25.0_r726/hifiasm/hifiasm"


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

def run(path_hifi_read, path_save, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):    

    cmd = f"""{PATH_HIFIASM} -o {path_save} --n-hap 1 -l0 --telo-m CCCTAA -t{n_jobs}"""
    
    for param_name, param_value in additional_params.items():
        if len(param_name) == 1:
            cmd += f" -{param_name} {param_value}"
        else:
            cmd += f" --{param_name} {param_value}"
    cmd += f" {path_hifi_read}"
        
    os.makedirs(os.path.dirname(path_save), exist_ok=True)
    os.makedirs(dir_tmp, exist_ok=True)
    os.makedirs(dir_shell, exist_ok=True)
    os.makedirs(dir_log, exist_ok=True)
    
    # print(cmd)
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
    
    list_contignames = [
        "h1tg000011l",
        "h1tg000012l",
        "h1tg000026l"
    ]
    
    for contig in list_contignames:
        
        path_hifi = f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/RawReads_Aligned_HiFi/from_round1/{contig}_elongate_1_end.R1.concise.hifi_readid.fastq"

        path_save = f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Assembly_Elongation/from_round1/{contig}_elongate_1_end.Hifiasm.HiFi/{contig}_elongate_1_end.hifi.haploid.l0.asm"

        assert os.path.exists(path_hifi), path_hifi

        curr_time = int(time.time())
        run(
            path_hifi,  
            path_save,
            f"{dir_workspace}/Hifiasm_Local_Assembly/tmp",
            f"{dir_workspace}/Hifiasm_Local_Assembly/sh",
            f"{dir_workspace}/Hifiasm_Local_Assembly/log",
            f"Hifiasm_Local_Assembly.HiFi.{contig}_elongate_1_end.l0.{curr_time}",
            64,
            hostname = [],
            jobs_wait = [],
            # **{
            #     "ul" : path_ont_hq_all
            # }
        )