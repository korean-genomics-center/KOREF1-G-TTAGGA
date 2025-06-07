import os, time
import subprocess

PATH_NTLINK = "/BiO/Access/yoonsung/anaconda3/envs/LongStitch/bin/ntLink"


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

def run(path_draft, path_reads, path_output_prefix, k_ntlink, window_ntlink, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), rounds = 1, gapfilling = False, **additional_params):
    dir_save = os.path.dirname(path_output_prefix)
    output_prefix = os.path.basename(path_output_prefix)
    os.makedirs(os.path.dirname(path_output_prefix), exist_ok=True)

    draft_dst = make_symlink(path_draft, dir_save, ch_extension="fa")
    read_dst = make_symlink(path_reads, dir_save)

    if rounds == 1:
        cmd = f"""{PATH_NTLINK} -C {dir_save} scaffold {'gap-fill' if gapfilling else ''} target={draft_dst} reads={read_dst} prefix={output_prefix} t={n_jobs} k={k_ntlink} w={window_ntlink}"""
    
    else:
        cmd = f"""{PATH_NTLINK}_rounds -C {dir_save} {'run_rounds_gaps' if gapfilling else 'run_rounds'} target={draft_dst} reads={read_dst} t={n_jobs} k={k_ntlink} w={window_ntlink} rounds={rounds}"""
    
    for param_name, param_value in additional_params.items():
        if len(param_name) == 1:
            cmd += f" {param_name}={param_value}"
        else:
            cmd += f" {param_name}={param_value}"
    
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
    
    dict_hap_to_side = {
        "1" : "p",
        "2" : "m"
    }

    for hap in ["1", "2"]:
        path_draft = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/Telomere_Detection/Fasta_Split_by_Telomere/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg.no_telo.telo_single_masked.fasta"
        path_reads = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/ONT/Parental_Split/ONT_UL.Adapt_Filt.Edge_Trim_Front_1500_Tail_1500.Filter_L50kb_QV20_Q30Perc80.merged.triobin_split.triobin_{dict_hap_to_side[hap]}.fq.gz"
        path_output_prefix = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.rounds5_gapfilling/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.no_telo.telo_single_masked.ntLink.ONT_triobin_{dict_hap_to_side[hap]}.L50kb_QV20_Q30Perc80.rounds5_gapfilling"
        k_ntlink = 32
        window_ntlink = 100
        rounds = 5
        gapfilling = True

        modified_option = {
        }

        run(
            path_draft, 
            path_reads, 
            path_output_prefix, 
            k_ntlink, 
            window_ntlink,
            f"{dir_workspace}/ntLink_Scaffolding/tmp",
            f"{dir_workspace}/ntLink_Scaffolding/sh",
            f"{dir_workspace}/ntLink_Scaffolding/log",
            f"ntLink.Hifiasm_Assembly.KOREF_S1.ul_L100kb_QV15.Trio_Specific.hap{hap}.with_ONT_L50kb_QV20_Q30Perc80_Triobin_{dict_hap_to_side[hap]}.rounds5_gapfilling.{curr_time}",
            24,
            [],
            [],
            rounds = rounds, 
            gapfilling = gapfilling,
            **modified_option
        )
#%%