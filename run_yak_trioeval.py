import os, time
import subprocess

PATH_YAK = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/yak/yak/yak"


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
    # subprocess.run(f"qsub {path_shell}", shell = True)
    print(path_shell)

def run(path_assembly, path_bin_pat, path_bin_mat, path_save_out, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):
    cmd_filtervariant = f"""{PATH_YAK} trioeval -e -t{n_jobs} {path_bin_pat} {path_bin_mat} {path_assembly} > {path_save_out}"""
    
    for param_name, param_value in additional_params.items():
        cmd_filtervariant += f" --{param_name} {param_value}"
    
    os.makedirs(os.path.dirname(path_save_out), exist_ok=True)
    os.makedirs(dir_tmp, exist_ok=True)
    os.makedirs(dir_shell, exist_ok=True)
    os.makedirs(dir_log, exist_ok=True)
    
    write_shell_script_and_run(
        f"{dir_shell}/{jobname}.sh",
        get_context_for_qsub_submission(
            cmd_filtervariant,
            jobname,
            f"{dir_log}/{jobname}.stderr",
            f"{dir_log}/{jobname}.stdout",
            n_jobs,
            hostname=hostname,
            hold_jids=jobs_wait
        )
    )



if __name__ == "__main2__":
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis"
    
    
    curr_time = int(time.time())
    
    for hap in ["1", "2"]:
        path_assembly = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg.fasta"
        paternal_count = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/PATERNAL/Yak_Binning/PATERNAL_Trimmed_F15_T5.noNbase.merged.yak"
        maternal_count = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/MATERNAL/Yak_Binning/MATERNAL_Trimmed_F15_T5.noNbase.merged.yak"
        
        path_save_stdout = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TrioEval/KKOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg.trioeval_err_stdout"
        
        curr_time = int(time.time())
        run(
            path_assembly,
            paternal_count,
            maternal_count,
            path_save_stdout,
            f"{dir_workspace}/Evaluate_Assembly_Triobin/tmp",
            f"{dir_workspace}/Evaluate_Assembly_Triobin/sh",
            f"{dir_workspace}/Evaluate_Assembly_Triobin/log",
            f"Evaluate_Assembly_Triobin.KOREF.ul_L100kb_QV15.Trio_Specific.hap{hap}.{curr_time}",
            16,
            [],
            []
        )


if __name__ == "__main__":
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis"
    
    dict_hap_to_triobin = {
        "1":'p',
        "2":'m'
    }

    curr_time = int(time.time())
    
    for hap in ["1", "2"]:
        path_assembly = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/ONT_UL.Filter_L50kb_QV20_Q30Perc80.with_0/Hap{hap}/assembly.fasta"
        paternal_count = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/Binning_Compare/Fasta_from_Meryl_Print/Yak_Binning/PATERNAL_Trimmed_F15_T5.noNbase.merged.k31.only.kmer_count.b0_k31.yak"
        maternal_count = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/Binning_Compare/Fasta_from_Meryl_Print/Yak_Binning/MATERNAL_Trimmed_F15_T5.noNbase.merged.k31.only.kmer_count.b0_k31.yak"
        
        path_save_stdout = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/ONT_UL.Filter_L50kb_QV20_Q30Perc80.with_0/Hap{hap}/assembly.parental_only_kmer.trioeval_err_stdout"
        
        curr_time = int(time.time())
        run(
            path_assembly,
            paternal_count,
            maternal_count,
            path_save_stdout,
            f"{dir_workspace}/Evaluate_Assembly_Triobin/tmp",
            f"{dir_workspace}/Evaluate_Assembly_Triobin/sh",
            f"{dir_workspace}/Evaluate_Assembly_Triobin/log",
            f"Evaluate_Assembly_Triobin.KOREF.ul_L100kb_QV15.Trio_Specific.hap{hap}.parental_only_kmer.{curr_time}",
            16,
            [],
            []
        )