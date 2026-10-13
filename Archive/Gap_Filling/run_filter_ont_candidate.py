import pandas as pd
import os, time
import subprocess


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

def run(path_fastq_in, path_fastq_out, path_paf, align_ratio, identity_ratio, aligned_region, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list()):    

    cmd = f"""python /BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Resources/Scripts/Gap_Filling/filter_ont_reads_for_gapclosing.py --fastq_in {path_fastq_in} --fastq_out {path_fastq_out} --paf {path_paf} --align_ratio {align_ratio} --identity_ratio {identity_ratio} --aligned_region {aligned_region}"""
        
    os.makedirs(os.path.dirname(path_fastq_out), exist_ok=True)
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
      
    dict_hap_to_triobin = {
        1 : 'p',
        2 : 'm'
    }
    
    for hap, triobin in dict_hap_to_triobin.items():
        
        path_fastq_in = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/ONT/Parental_Split/ONT_UL.Adapt_Filt.Edge_Trim_Front_1500_Tail_1500.Filter_L50kb_QV20_Q30Perc80.merged.triobin_split.triobin_{triobin}.fq.gz"

        path_fastq_out = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/Round2.Unplaced_Scf_asm5X_idy5_match1000.ONT_UL_L50kb_QV20_Q30Perc80_parental/Pre_Alignment.Filter_Seen_Sequence/Hap{hap}/n_pseudo_telorepeat5/Candidate_Sequence.ONT_UL.Adapt_Filt.Edge_Trim_Front_1500_Tail_1500.Filter_L50kb_QV20_Q30Perc80.merged.triobin_split.triobin_{triobin}.filterout_align_0.95_ident_0.8_region_1Mb.fq.gz"
        
        path_paf = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/Round2.Unplaced_Scf_asm5X_idy5_match1000.ONT_UL_L50kb_QV20_Q30Perc80_parental/Pre_Alignment.Filter_Seen_Sequence/Hap{hap}/n_pseudo_telorepeat5/Minimap2_Align.Ref_CuratedDraftASM.Trio_Specific.hap{hap}.n_pseudo_telorepeat5.GapFilling.Unplaced_Scf.SeqSplit.Qry_ONT_UL.Filter_L50kb_QV20_Q30Perc80.triobin_{triobin}.mapont.paf"

        align_ratio = 0.95
        identity_ratio = 0.8
        aligned_region = 1000000

        curr_time = int(time.time())
        run(
            path_fastq_in, 
            path_fastq_out, 
            path_paf, 
            align_ratio, 
            identity_ratio, 
            aligned_region, 
            f"{dir_workspace}/Filter_Gapclosing_Candidate_Sequence/tmp",
            f"{dir_workspace}/Filter_Gapclosing_Candidate_Sequence/sh",
            f"{dir_workspace}/Filter_Gapclosing_Candidate_Sequence/log",
            f"Filter_Gapclosing_Candidate_Sequence.ONT_UL_L50kb_QV20_Q30Perc80_parental.triobin_{triobin}.{curr_time}",
            1,
            hostname = [],
            jobs_wait = []
        )