#%%
import pandas as pd
import os, time
import subprocess

PATH_SCRIPT = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Resources/Scripts/remove_duplicated_fastq_reads.py"


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

def run_dedup(path_fastq, path_save, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):
    cmd_filtervariant = f"""python {PATH_SCRIPT} --input {path_fastq} --out {path_save}"""
    
    for param_name, param_value in additional_params.items():
        cmd_filtervariant += f" --{param_name} {param_value}"
    
    os.makedirs(os.path.dirname(path_save), exist_ok=True)
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

if __name__ == "__main__":
    list_paths = [
        "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/Align.Winnowmap2/Region_filtered/BAM2FastQ.Ref.CHM13.QRY.ONT_UL.Filter_L50kb_QV20_Q30Perc80.triobin_split.triobin_m.chr15_98753195.fastq",
        "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/Align.Winnowmap2/Region_filtered/BAM2FastQ.Ref.CHM13.QRY.ONT_UL.Filter_L50kb_QV20_Q30Perc80.triobin_split.triobin_p.chr13_1_1000000.fastq",
        "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/Align.Winnowmap2/Region_filtered/BAM2FastQ.Ref.CHM13.QRY.ONT_UL.Filter_L50kb_QV20_Q30Perc80.triobin_split.triobin_p.chr13_1_20000000.fastq",
        "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/Align.Winnowmap2/Region_filtered/BAM2FastQ.Ref.CHM13.QRY.ONT_UL.Filter_L50kb_QV20_Q30Perc80.triobin_split.triobin_p.chrY_20000000.fastq",
        "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/Align.Winnowmap2/Region_filtered/BAM2FastQ.Ref.CHM13.QRY.PacBio_HiFi_ccs_converted.Filter_L1000_QV20.triobin_split.triobin_m.chr15_98753195.fastq",
        "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/Align.Winnowmap2/Region_filtered/BAM2FastQ.Ref.CHM13.QRY.PacBio_HiFi_ccs_converted.Filter_L1000_QV20.triobin_split.triobin_p.chr13_1_1000000.fastq",
        "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/Align.Winnowmap2/Region_filtered/BAM2FastQ.Ref.CHM13.QRY.PacBio_HiFi_ccs_converted.Filter_L1000_QV20.triobin_split.triobin_p.chr13_1_20000000.fastq",
        "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/Align.Winnowmap2/Region_filtered/BAM2FastQ.Ref.CHM13.QRY.PacBio_HiFi_ccs_converted.Filter_L1000_QV20.triobin_split.triobin_p.chrY_20000000.fastq"
    ]

    dir_save = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/Align.Winnowmap2/Region_filtered/Deduplication"
    
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis/Dedup_FQ"
    
    for path_input in list_paths:
        filename = os.path.basename(path_input)
        filename_save = filename.replace(".fastq", ".dedup.fastq")
        filename_save_unsafe = filename.replace(".fastq", ".unsafe_dup.fastq")
        
        path_save = os.path.join(dir_save, filename_save)
        path_save_unsafe = os.path.join(dir_save, filename_save_unsafe)
        
        run_dedup(
                path_input,
                path_save,
                f"{dir_workspace}/Dedup_Fastq/tmp",
                f"{dir_workspace}/Dedup_Fastq/sh",
                f"{dir_workspace}/Dedup_Fastq/log",
                f"Dedup_Fastq.Local_Assembly.{filename.replace('.fastq', '')}.{int(time.time())}",
                1,
                [],
                unsafe = path_save_unsafe
        )