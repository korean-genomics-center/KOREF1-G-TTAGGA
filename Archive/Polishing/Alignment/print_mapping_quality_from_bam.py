import pandas as pd
import os, time
import subprocess

PATH_SAMTOOLS = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/Samtools/samtools-1.17/samtools-1.17/bin/samtools"


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

def run(path_input_bam, path_output, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list()):    

    cmd = """{PATH_SAMTOOLS} view -@ {n_jobs} {path_input_bam} | awk '{read_name = $1;mapq = $5;seq = $10;identity_str = "NA";nm_tag = -1; for (i=12; i<=NF; ++i) {if ($i ~ /^NM:i:/) {nm_tag = substr($i, 6);break;}} seq_len = length(seq); if (seq_len > 0 && nm_tag != -1) {identity = (seq_len - nm_tag) / seq_len * 100; identity_str = sprintf("%.4f", identity);} print read_name, mapq, seq_len, identity_str;}' > {path_output}""".replace("{PATH_SAMTOOLS}", PATH_SAMTOOLS).replace("{n_jobs}", str(n_jobs)).replace("{path_input_bam}", path_input_bam).replace("{path_output}", path_output)
        
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

if __name__ == "__main2__":
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis"
    
    for hap in [1, 2]:
        for parental in ['m', 'p']:
        # for parental in ['a', '0', 'm', 'p']:
                
            # path_input_bam = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Read_Alignment/Align.Winnowmap2/Winnowmap_Align.Ref.KOREF_ASM.hap{hap}.QRY.PacBio_HiFi_ccs_converted.Filter_L1000_QV20.triobin_split.triobin_{parental}.bam"
            # path_output_bam = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Read_Alignment/Align.Winnowmap2/Primary_Alignment/Winnowmap_Align.Primary.Ref.KOREF_ASM.hap{hap}.QRY.PacBio_HiFi_ccs_converted.Filter_L1000_QV20.triobin_split.triobin_{parental}.bam"
            
            path_input_bam = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Read_Alignment/Align.Winnowmap2/Align_Gapclosed_ASM/Primary_Alignment/Winnowmap_Align.Primary.Ref.KOREF_ASM.hap{hap}.QRY.ONT_UL.Filter_L10kb_QV15.triobin_split.triobin_{parental}.bam"
            path_output = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Read_Alignment/Align.Winnowmap2/Align_Gapclosed_ASM/Primary_Alignment/Winnowmap_Align.Primary.Ref.KOREF_ASM.hap{hap}.QRY.ONT_UL.Filter_L10kb_QV15.triobin_split.triobin_{parental}.mapq_identity"
            
            curr_time = int(time.time())
            run(
                path_input_bam, 
                path_output, 
                f"{dir_workspace}/Polishing/Check_MapQ_Identity/tmp", 
                f"{dir_workspace}/Polishing/Check_MapQ_Identity/sh", 
                f"{dir_workspace}/Polishing/Check_MapQ_Identity/log", 
                f"Check_MapQ_Identity.ONT_UL.Filter_L10kb_QV15.triobin_split.REF.KOREF_ASM.hap{hap}.triobin_{parental}.{curr_time}",
                8
            )
            
if __name__ == "__main__":
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis"
    
    for hap in [1, 2]:
        # for parental in ['m', 'p']:
        for parental in ['a', '0', 'm', 'p']:
            
            path_input_bam = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Read_Alignment/Align.Winnowmap2/Align_Gapclosed_ASM/Primary_Alignment/Winnowmap_Align.Primary.Ref.KOREF_ASM.hap{hap}.QRY.PacBio_HiFi_ccs_converted.Filter_L1000_QV20.triobin_split.triobin_{parental}.bam"
            path_output = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Read_Alignment/Align.Winnowmap2/Align_Gapclosed_ASM/Primary_Alignment/Winnowmap_Align.Primary.Ref.KOREF_ASM.hap{hap}.QRY.PacBio_HiFi_ccs_converted.Filter_L1000_QV20.triobin_split.triobin_{parental}.mapq_identity"
            
            curr_time = int(time.time())
            run(
                path_input_bam, 
                path_output, 
                f"{dir_workspace}/Polishing/Check_MapQ_Identity/tmp", 
                f"{dir_workspace}/Polishing/Check_MapQ_Identity/sh", 
                f"{dir_workspace}/Polishing/Check_MapQ_Identity/log", 
                f"Check_MapQ_Identity.PacBio_HiFi_ccs_converted.Filter_L1000_QV20.triobin_split.REF.KOREF_ASM.hap{hap}.triobin_{parental}.{curr_time}",
                8
            )