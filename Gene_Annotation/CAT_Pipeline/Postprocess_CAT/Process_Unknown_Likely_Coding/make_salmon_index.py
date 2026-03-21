#%%
import subprocess, os, time

PATH_SALMON = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/Salmon/Salmon_v1.10.0/salmon-latest_linux_x86_64/bin/salmon"

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

def run(path_transcript, path_genome, path_output_prefix, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):
    dir_output = os.path.dirname(path_output_prefix)

    cmd = f"""
grep "^>" {path_genome} | cut -d " " -f 1 > {path_output_prefix}.decoys.txt
sed -i.bak -e 's/>//g' {path_output_prefix}.decoys.txt
cat {path_transcript} {path_genome} > {path_output_prefix}.gentrome.fa
{PATH_SALMON} index -t {path_output_prefix}.gentrome.fa -d {path_output_prefix}.decoys.txt -p {n_jobs} -i {path_output_prefix}.salmon_index
"""
    
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
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis"

    curr_time = int(time.time())

    for hap in [1, 2]:
        path_transcript = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/CDS_Fasta/Hap{hap}.Blood_Cellline.CGP_PB/KOREF1_ASM_Hap{hap}.whole_gff.transcript.fa"
        path_genome = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/CDS_Fasta/Hap{hap}.Blood_Cellline.CGP_PB/Salmon/Index_Salmon/KOREF1_ASM_Hap{hap}.decoy.fa"

        path_output_prefix = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/CDS_Fasta/Hap{hap}.Blood_Cellline.CGP_PB/Salmon/Index_Salmon/KOREF1_ASM_Hap{hap}.whole_gff"
        
        run(
            path_transcript, path_genome, path_output_prefix,
            f"{dir_workspace}/CAT_Annotation/Run_Salmon_Index/tmp",
            f"{dir_workspace}/CAT_Annotation/Run_Salmon_Index/sh",
            f"{dir_workspace}/CAT_Annotation/Run_Salmon_Index/log",
            f"Run_Salmon_Index.hap{hap}.{curr_time}",
            64,
            [],
            []
        )