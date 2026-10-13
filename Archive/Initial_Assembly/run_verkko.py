import os, time
import subprocess

PATH_VERKKO = "/BiO/Access/yoonsung/anaconda3/envs/verkko/bin/verkko"
PATH_GRAPHALIGNER = "/BiO/Access/yoonsung/anaconda3/envs/verkko_dependencies_new/bin/GraphAligner"
PATH_MASHMAP = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/Mashmap/Mashmap_v3.1.3/MashMap/build/bin/mashmap"
PATH_WINNOWMAP = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/winnowmap/winnowmap.v2.03/Winnowmap/bin/winnowmap"

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

#--graphaligner {PATH_GRAPHALIGNER} --mashmap {PATH_MASHMAP} --winnowmap {PATH_WINNOWMAP} 
def run(path_hifi_read, path_ul_read, paternal_count, maternal_count, dir_output, dir_tmp, dir_shell, dir_log, jobname, n_jobs, memory, hostname = list(), jobs_wait = list(), **additional_params):
    cmd = f"""{PATH_VERKKO} --mashmap {PATH_MASHMAP} --winnowmap {PATH_WINNOWMAP} -d {dir_output} --hifi {path_hifi_read} --nano {path_ul_read} --hap-kmers {paternal_count} {maternal_count} trio --local-memory {memory} --local-cpus {n_jobs}"""
    
    for param_name, param_value in additional_params.items():
        if len(param_name) == 1:
            cmd += f" -{param_name} {param_value}"
        else:
            cmd += f" --{param_name} {param_value}"
    
    os.makedirs(dir_output, exist_ok=True)
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
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis/Verkko_Assembly"
    
    
    curr_time = int(time.time())
    
    path_hifi_read = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/PacBio/Merged/PacBio_HiFi_ccs_converted.filt.Filter_L1000_QV20.fastq.gz"
    path_ul_read = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/ONT/Merged/ONT_UL.Adapt_Filt.Edge_Trim_Front_1500_Tail_1500.Filter_L100kb_QV15.merged.fastq"
    paternal_count = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/PacBio/Hapmer_Diff/Hapmer_diff_k30.compress/PATERNAL_Trimmed_F15_T5.noNbase.merged.R12.k30.compress.hapmer.meryl"
    maternal_count = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/PacBio/Hapmer_Diff/Hapmer_diff_k30.compress/MATERNAL_Trimmed_F15_T5.noNbase.merged.R12.k30.compress.hapmer.meryl"
    dir_output = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Verkko/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.default_option"

    n_jobs = 128
    memory = 1000
    
    curr_time = int(time.time())
    
    run(
        path_hifi_read, 
        path_ul_read, 
        paternal_count, 
        maternal_count, 
        dir_output,
        f"{dir_workspace}/KOREF_S1_TrioBinning/tmp",
        f"{dir_workspace}/KOREF_S1_TrioBinning/sh",
        f"{dir_workspace}/KOREF_S1_TrioBinning/log",
        f"Verkko_Assembly.KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.default_option.{curr_time}",
        n_jobs,
        memory,
        [],
        [],
        # **moditifed_option
    )
#%%