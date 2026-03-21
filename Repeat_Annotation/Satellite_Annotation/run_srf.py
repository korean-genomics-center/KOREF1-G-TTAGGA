import os, time, math
import subprocess
from copy import deepcopy

PATH_KMC = "/BiO/Research/Project2/Endangered_Jindo_Dog/Resources/Tools/KMC/KMC_v3.2.4/bin/kmc"
PATH_KMC_DUMP = "/BiO/Research/Project2/Endangered_Jindo_Dog/Resources/Tools/KMC/KMC_v3.2.4/bin/kmc_dump"
PATH_SRF = "/BiO/Research/Project2/Endangered_Jindo_Dog/Resources/Tools/SRF/SRF_commit_e54ca8c/srf/srf"


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

def run(path_asm, path_output_prefix, kmer_size, min_k_count, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):
    path_tmp_dir = f"{path_output_prefix}.kmc_tmp"
    cmd = f"""
{PATH_KMC} -fa -k{kmer_size} -t{n_jobs} -ci{min_k_count} -cs100000 {path_asm} {path_output_prefix}.kmer_count.kmc {path_tmp_dir};
{PATH_KMC_DUMP} {path_output_prefix}.kmer_count.kmc {path_output_prefix}.kmer_count.txt;
{PATH_SRF} -p repeat {path_output_prefix}.kmer_count.txt > {path_output_prefix}.srf_output.fa;
"""
           
    os.makedirs(os.path.dirname(path_output_prefix), exist_ok=True)
    os.makedirs(path_tmp_dir, exist_ok=True)
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
    dir_workspace = "/BiO/Research/Project2/Endangered_Jindo_Dog/Analysis"

    curr_time = int(time.time())

    modified_option = {
    }
    
    kmer_size = 101
    min_k_count = 10
    
    for hap in [1, 2]:
        path_asm = f"/BiO/Research/Project2/Endangered_Jindo_Dog/Results/Repeat_Annotation/Satellite_Repeat/run_srf/Hap{hap}/FinalASM.J495799_Child.Trio_Specific.hap{hap}.fasta"
        path_output_prefix = f"/BiO/Research/Project2/Endangered_Jindo_Dog/Results/Repeat_Annotation/Satellite_Repeat/run_srf/Hap{hap}/Satellite_repeat_detection.J495799_Child.Hap{hap}.k{kmer_size}.ci{min_k_count}"
        
        jobname = f"J495799_Child_Hap{hap}.k{kmer_size}.ci{min_k_count}"
        
        run(
            path_asm, path_output_prefix, kmer_size, min_k_count,
            f"{dir_workspace}/Repeat_Annotation/Satellite_Detection/tmp",
            f"{dir_workspace}/Repeat_Annotation/Satellite_Detection/sh",
            f"{dir_workspace}/Repeat_Annotation/Satellite_Detection/log",
            f"Satellite_Detection.{jobname}.{curr_time}",
            32,
            [],
            [],
            **modified_option
        )
        
    
    # path_asm = f"/BiO/Research/Project2/Endangered_Jindo_Dog/Results/Repeat_Annotation/Satellite_Repeat/run_srf/Mixed/FinalASM.J495799_Child.Trio_Specific.hap_Mixed.fasta"
    # path_output_prefix = f"/BiO/Research/Project2/Endangered_Jindo_Dog/Results/Repeat_Annotation/Satellite_Repeat/run_srf/Mixed/Satellite_repeat_detection.J495799_Child.Hap_Mixed.k{kmer_size}.ci{min_k_count}"
    
    # jobname = f"J495799_Child_Hap_Mixed.k{kmer_size}.ci{min_k_count}"
    
    # run(
    #     path_asm, path_output_prefix, kmer_size, min_k_count,
    #     f"{dir_workspace}/Repeat_Annotation/Satellite_Detection/tmp",
    #     f"{dir_workspace}/Repeat_Annotation/Satellite_Detection/sh",
    #     f"{dir_workspace}/Repeat_Annotation/Satellite_Detection/log",
    #     f"Satellite_Detection.{jobname}.{curr_time}",
    #     64,
    #     [],
    #     [],
    #     **modified_option
    # )
#%%