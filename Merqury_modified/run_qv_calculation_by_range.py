#%%
import subprocess, os, time

PATH_SCRIPT = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Resources/Scripts/Merqury_modified/per_seq_qv_range.sh"
DIR_MERQURY = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/merqury"
DIR_MERYL = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/meryl-1.4.1/bin"
DIR_SAMTOOLS = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/Samtools/samtools-1.17/samtools-1.17/bin"
DIR_BEDTOOLS = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/Bedtools"

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
        os.remove(path_dst)
    os.symlink(path_src, path_dst)

def run(path_asm, path_child_meryl, path_output_prefix, window_size, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):
    path_dir_save = os.path.dirname(path_output_prefix)
    output_prefix = os.path.basename(path_output_prefix)

    path_asm_link = os.path.join(path_dir_save, os.path.basename(path_asm))
    os.makedirs(path_dir_save, exist_ok=True)
    make_symlink(path_asm, path_asm_link)

    cmd = f"""cd {path_dir_save};{PATH_SCRIPT} {path_asm_link} {path_child_meryl} {output_prefix} {window_size}"""
    
    merqury_environ = os.environ.copy()
    merqury_environ["PATH"] = f"{DIR_MERYL}:{DIR_SAMTOOLS}:{DIR_BEDTOOLS}:{merqury_environ['PATH']}"
    merqury_environ["MERQURY"] = DIR_MERQURY
    merqury_environ["OMP_NUM_THREADS"] = f"{n_jobs}"
    
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
        ),
        env = merqury_environ
    )

if __name__ == "__main__":
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis"
    
    k = 21

    path_asm_format = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Reference_Fasta/Round8.SV_SNV_Polished_Fasta/PolishedASM.KOREF_S1.Trio_Specific.UpperCase.Polish_R8.SV_SNV.hap{hap}.fasta"
    
    path_child_meryl = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/PacBio/Meryl_Binning/PacBio_HiFi_ccs_converted.filt.Filter_L1000_QV20.k{k}.meryl"
    
    path_output_prefix_format = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Reference_Fasta/Round8.SV_SNV_Polished_Fasta/Merqury_QV_per_window/PolishedASM.KOREF_S1.Trio_Specific.UpperCase.Polish_R8.SV_SNV.hap{hap}.QV_per_window"

    window_size = 10000
    
    curr_time = int(time.time())

    for hap in [1, 2]:
        run(
            path_asm_format.format(hap = hap),
            path_child_meryl,
            path_output_prefix_format.format(hap = hap, k = k),
            window_size,
            f"{dir_workspace}/Merqury_Evaluation/tmp",
            f"{dir_workspace}/Merqury_Evaluation/sh",
            f"{dir_workspace}/Merqury_Evaluation/log",
            f"Evaluate_QV_per_window.PolishedASM.KOREF_S1.Polish_R8.SV_SNV.Hap{hap}.k{k}.{window_size}.{curr_time}",
            64,
            [],
            []
        )
#%%