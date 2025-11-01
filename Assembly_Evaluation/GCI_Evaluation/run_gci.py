import pandas as pd
import os, time
import subprocess

PATH_GCI = "/BiO/Access/yoonsung/anaconda3/envs/gci/bin/GCI"


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

def run(path_reference, path_hifi_winnowmap2_bam, path_hifi_minimap2_paf, path_ont_winnowmap2_bam, path_ont_minimap2_paf, path_output_prefix, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list()):    
    dir_out = os.path.dirname(path_output_prefix)
    output_prefix = os.path.basename(path_output_prefix)

    cmd = f"""{PATH_GCI} -r {path_reference} --hifi {path_hifi_winnowmap2_bam} {path_hifi_minimap2_paf} --nano {path_ont_winnowmap2_bam} {path_ont_minimap2_paf} -d {dir_out} -o {output_prefix} -t {n_jobs} -p -it pdf"""
        
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
    dict_hap_to_parental = {
        1:'p',
        2:'m'
    }
    
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis"

    for hap in [1, 2]:
        
        path_reference = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Final_Assembly/FinalASM.KOREF_S1.Trio_Specific.hap{hap}.Reordered.fasta"
        

        assert os.path.exists(path_fasta_qry)
        assert os.path.exists(path_fasta_ref)

        
        curr_time = int(time.time())
        run(
            path_reference,
            path_hifi_winnowmap2_bam,
            path_hifi_minimap2_paf,
            path_ont_winnowmap2_bam,
            path_ont_minimap2_paf,
            path_output_prefix,
            f"{dir_workspace}/GCI_Evaluation/Run_GCI/tmp",
            f"{dir_workspace}/GCI_Evaluation/Run_GCI/sh",
            f"{dir_workspace}/GCI_Evaluation/Run_GCI/log",
            f"Run_GCI.REF.KOREF_FinalASM.hap{hap}.{curr_time}",
            64,
            hostname = [],
            jobs_wait = []
        )