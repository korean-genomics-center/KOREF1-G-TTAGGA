#%%
import subprocess, os, time

PATH_SAMTOOLS = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/Samtools/samtools-1.17/samtools-1.17/bin/samtools"
PATH_BWA = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/bwa/bwa-0.7.19/bwa/bwa"
PATH_PAIRTOOLS = "/BiO/Access/yoonsung/anaconda3/envs/pairtools/bin/pairtools"


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

def run(path_asm, path_hic_1, path_hic_2, path_output_prefix, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):
    dir_save = os.path.dirname(path_output_prefix)
    
    path_copy_asm = os.path.join(dir_save, os.path.basename(path_asm))
    path_save_asm_size = os.path.join(dir_save, f"{os.path.basename(path_asm)}.genome")
    
    dir_output_tmp = f"{dir_save}/pairtools_tmp"
    
    cmd_preprocessing = f"""
cp {path_asm} {path_copy_asm};
{PATH_SAMTOOLS} faidx {path_copy_asm};
cut -f1,2 {path_copy_asm}.fai > {path_save_asm_size};
{PATH_BWA} index {path_copy_asm};
"""

    cmd_process_hic = f"""
{PATH_BWA} mem -5SP -T0 -t{n_jobs} {path_copy_asm} {path_hic_1} {path_hic_2} | {PATH_SAMTOOLS} view -bS -@ {min(n_jobs, 8)} -o {path_output_prefix}.aligned.bam;
{PATH_PAIRTOOLS} parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in {n_jobs} --nproc-out {n_jobs} --chroms-path {path_save_asm_size} -o {path_output_prefix}.pairsam.gz {path_output_prefix}.aligned.bam
{PATH_PAIRTOOLS} sort --nproc {n_jobs} --tmpdir={dir_output_tmp} {path_output_prefix}.pairsam.gz -o {path_output_prefix}.sorted.pairsam.gz;
{PATH_PAIRTOOLS} dedup --nproc-in {n_jobs} --nproc-out {n_jobs} --mark-dups --output-stats {path_output_prefix}.sorted.dedup.pairsam.stats --output {path_output_prefix}.sorted.dedup.pairsam.gz {path_output_prefix}.sorted.pairsam.gz;
{PATH_PAIRTOOLS} split --nproc-in {n_jobs} --nproc-out {n_jobs} --output-pairs {path_output_prefix}.split.pairs.gz --output-sam {path_output_prefix}.split.bam {path_output_prefix}.sorted.dedup.pairsam.gz;
{PATH_SAMTOOLS} sort -@{n_jobs} -o {path_output_prefix}.split.sorted.bam {path_output_prefix}.split.bam;
{PATH_SAMTOOLS} index -@{n_jobs} {path_output_prefix}.split.sorted.bam;
"""

    cmd = ""
    if not os.path.exists(path_copy_asm):
        cmd += cmd_preprocessing
    cmd += cmd_process_hic
    
    os.makedirs(dir_save, exist_ok=True)
    os.makedirs(dir_output_tmp, exist_ok=True)
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
    
    for hap in [1, 2]:
        path_asm = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Reference_Fasta/Round8.SV_SNV_Polished_Fasta/PolishedASM.KOREF_S1.Trio_Specific.UpperCase.Polish_R8.SV_SNV.hap{hap}.fasta"
        
        path_hic_1 = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/HiC/Trim_FQ/KOREF-S1_PBMC_Hi-C_NovaSeq_20190425_DpnII-HinfI_L01_R1.dedup.trim_F15_T0.fq.gz"
        path_hic_2 = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/HiC/Trim_FQ/KOREF-S1_PBMC_Hi-C_NovaSeq_20190425_DpnII-HinfI_L01_R2.dedup.trim_F15_T0.fq.gz"
        
        path_output_prefix = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/HiC_Evaluation/HiC_Alignment/Hap{hap}/HiC_Processing.Ref.KOREF_S1.PolishedASM_Hap{hap}"
        
        curr_time = int(time.time())
        run(
            path_asm,
            path_hic_1,
            path_hic_2,
            path_output_prefix,
            f"{dir_workspace}/Evaluate_Assembly/HiC_Alignment/tmp",
            f"{dir_workspace}/Evaluate_Assembly/HiC_Alignment/sh",
            f"{dir_workspace}/Evaluate_Assembly/HiC_Alignment/log",
            f"HiC_Alignment.KOREF_S1.PolishedASM_Hap.Hap{hap}.{curr_time}",
            128,
            [],
            []
        )
#%%