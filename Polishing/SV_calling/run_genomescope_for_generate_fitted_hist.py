import pandas as pd
import os, time
import subprocess

PATH_MERYL = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/meryl-1.4.1/bin/meryl"
PATH_RSCRIPT = "/BiO/Access/yoonsung/anaconda3/envs/genomescope2/bin/Rscript"
PATH_GENOMESCOPE2 = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/GenomeScope2/GenomeScope_v2.0.1/genomescope2.0/genomescope.R"


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

def run(dir_meryldb, kmer_size, path_output_prefix, dir_tmp, dir_shell, dir_log, jobname, n_jobs, run_meryl = True, hostname = list(), jobs_wait = list(), genomescope2_params = dict()):  
    dir_save_gt1_meryldb = f"{path_output_prefix}.readmer.gt1.meryl"
    path_save_merylhist = f"{path_output_prefix}.readmer.hist"
    dir_save_genomescope = f"{path_output_prefix}.readmer.genomescope_out"

    if run_meryl:
        cmd = f"""{PATH_MERYL} greater-than 1 {dir_meryldb} output {dir_save_gt1_meryldb}; {PATH_MERYL} histogram {dir_meryldb} > {path_save_merylhist}; """
    else:
        cmd = ""

    genomescope2_param_line = ''
    for key, val in genomescope2_params.items():
        genomescope2_param_line += f" {key} {val}"

    cmd += f"{PATH_RSCRIPT} {PATH_GENOMESCOPE2} -i {path_save_merylhist} -k {kmer_size} -o {dir_save_genomescope} -p 2 --fitted_hist {genomescope2_param_line}"""
        
    os.makedirs(os.path.dirname(path_output_prefix), exist_ok=True)
    os.makedirs(dir_save_genomescope, exist_ok=True)
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
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis"
    
    kmer_size = 21

    dir_meryldb = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Read_Alignment/MerylDB/PacBio_HiFi_ccs_converted.filt.Filter_L1000_QV20.k{kmer_size}.meryl"
    path_output_prefix = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Read_Alignment/MerylDB/GenomeScope2/PacBio_HiFi_ccs_converted.filt.Filter_L1000_QV20.k{kmer_size}.l40"
    
    curr_time = int(time.time())
        
    run(
        dir_meryldb,
        kmer_size,
        path_output_prefix,
        f"{dir_workspace}/Polishing/GenomeScope_Seqmer/tmp",
        f"{dir_workspace}/Polishing/GenomeScope_Seqmer/sh",
        f"{dir_workspace}/Polishing/GenomeScope_Seqmer/log",
        f"GenomeScope_Seqmer.PacBio_HiFi_ccs_converted.filt.Filter_L1000_QV20.kmer{kmer_size}.{curr_time}",
        32,
        run_meryl = False,
        hostname = [],
        jobs_wait = [],
        genomescope2_params = {
            "-l" : "40"
        }
    )