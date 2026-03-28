import os, time
import subprocess

PATH_AGAT_KEEP = "/BiO/Access/yoonsung/anaconda3/envs/liftoff/bin/agat_sp_filter_feature_from_keep_list.pl"


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

def run(path_gff_in, path_keep_list, path_output, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list()):    

    cmd = f"""
{PATH_AGAT_KEEP} --gff {path_gff_in} --keep_list {path_keep_list} -o {path_output}

"""
        
    os.makedirs(os.path.dirname(path_output), exist_ok=True)
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
        for alert_level in ["filtout", "warn"]:
            path_gff_in = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Run_CAT/Hap{hap}.Blood_Cellline.CGP_PB/consensus_gene_set/KOREF1_ASM_Hap{hap}.gff3"
            path_keep_list = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/Process_Unknown/Hap{hap}.Blood_Cellline.CGP_PB/Evaluate.KOREF1_ASM_Hap{hap}.unknown_transcript.TPM_cutoff_1.transcrpts.{alert_level}.txt"
            path_output = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/Filter_QCfailed_Unknown/KOREF1_ASM_Hap{hap}.filteredout_unknown_transcript.TPM_cutoff_1.{alert_level}.gff3"

            assert os.path.exists(path_gff_in)
            assert os.path.exists(path_keep_list)
            run(
                path_gff_in, path_keep_list, path_output, 
                f"{dir_workspace}/CAT_Annotation/Filter_out_QCfailed/tmp", 
                f"{dir_workspace}/CAT_Annotation/Filter_out_QCfailed/sh", 
                f"{dir_workspace}/CAT_Annotation/Filter_out_QCfailed/log", 
                f"Get_QCfailed.KOREF1_FinalASM.CAT.hap{hap}.{curr_time}",
                1,
                jobs_wait=[]
            )