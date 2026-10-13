import pandas as pd
import os, time
import subprocess

PATH_HIFIASM = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/hifiasm_v0.25.0-r726/hifiasm/hifiasm"


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

def run_fastp_report(path_hifi_read, path_ul_read, paternal_count, maternal_count, path_save, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):
    cmd = f"""{PATH_HIFIASM} -o {path_save} --dual-scaf --telo-m CCCTAA -1 {paternal_count} -2 {maternal_count} --ul {path_ul_read} -t{n_jobs}"""
    
    for param_name, param_value in additional_params.items():
        if len(param_name) == 1:
            cmd += f" -{param_name} {param_value}"
        else:
            cmd += f" --{param_name} {param_value}"
    cmd += f" {path_hifi_read}"
    
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
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis/Hifiasm_Assembly"
    
    
    curr_time = int(time.time())
    
    path_hifi_read = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/PacBio/Merged/PacBio_HiFi_ccs_converted.filt.Filter_L1000_QV20.fastq.gz"
    path_ul_read = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/ONT/Merged/ONT_UL.Adapt_Filt.Edge_Trim_Front_1500_Tail_1500.Filter_L100kb_QV15.merged.fastq"
    paternal_count = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/Binning_Compare/Fasta_from_Meryl_Print/Yak_Binning/PATERNAL_Trimmed_F15_T5.noNbase.merged.k31.hapmer.kmer_count.b0_k31.yak"
    maternal_count = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/Binning_Compare/Fasta_from_Meryl_Print/Yak_Binning/MATERNAL_Trimmed_F15_T5.noNbase.merged.k31.hapmer.kmer_count.b0_k31.yak"
    path_save = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Hapmer.default_option/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Hapmer.asm"
    
    curr_time = int(time.time())
    
    # moditifed_option = {
    #     "z" : 20,
    #     "k" : 63,
    #     "D" : 7,
    #     "N" : 200,
    #     "r" : 5
    # }
    # moditifed_option = {
    #     "t-occ" : 5,
    #     'c' : 1,
    #     'd' : 10
    # }
    # moditifed_option = {
    #     "h1" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/HiC/Trim_FQ/KOREF-S1_PBMC_Hi-C_HiSeqXten_20150518_DpnII-HinfI_L01_R1.dedup.trim_F15_T0.fq.gz,/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/HiC/Trim_FQ/KOREF-S1_PBMC_Hi-C_HiSeqXten_20150518_MboI_L01_R1.dedup.trim_F15_T0.fq.gz,/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/HiC/Trim_FQ/KOREF-S1_PBMC_Hi-C_NovaSeq_20190425_DpnII-HinfI_L01_R1.dedup.trim_F15_T0.fq.gz",
    #     "h2" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/HiC/Trim_FQ/KOREF-S1_PBMC_Hi-C_HiSeqXten_20150518_DpnII-HinfI_L01_R2.dedup.trim_F15_T0.fq.gz,/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/HiC/Trim_FQ/KOREF-S1_PBMC_Hi-C_HiSeqXten_20150518_MboI_L01_R2.dedup.trim_F15_T0.fq.gz,/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/HiC/Trim_FQ/KOREF-S1_PBMC_Hi-C_NovaSeq_20190425_DpnII-HinfI_L01_R2.dedup.trim_F15_T0.fq.gz"
    # }
    
    run_fastp_report(
        path_hifi_read, 
        path_ul_read, 
        paternal_count, 
        maternal_count, 
        path_save,
        f"{dir_workspace}/KOREF_S1_TrioBinning/tmp",
        f"{dir_workspace}/KOREF_S1_TrioBinning/sh",
        f"{dir_workspace}/KOREF_S1_TrioBinning/log",
        f"Hifiasm_Assembly.KOREF_S1.ul_L100kb_QV15.Trio_Hapmer.default_option.{curr_time}",
        128,
        [],
        [],
        # **moditifed_option
    )
#%%