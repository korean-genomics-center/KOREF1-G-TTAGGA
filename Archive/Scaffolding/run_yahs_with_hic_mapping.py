import os, time
import subprocess

PATH_YAHS = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/YaHS/yahs_v1.2.2/yahs"
PATH_SAMTOOLS = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/Samtools/samtools-1.17/samtools-1.17/samtools"
PATH_BWAMEM = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/bwa/bwa-0.7.19/bwa/bwa"
PATH_PICARD = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/picard/picard-3.4.0/picard.jar"
PATH_JAVA = "/usr/bin/java"


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
    return path_dst

def run(path_draft, path_hic_read1, path_hic_read2, path_output_prefix, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), run_mapping = True, **additional_params):
    dir_save = os.path.dirname(path_output_prefix)
    os.makedirs(os.path.dirname(path_output_prefix), exist_ok=True)

    path_draft_lnk = make_symlink(path_draft, dir_save)

    cmd = ""
    if run_mapping:
        cmd += f"""
{PATH_SAMTOOLS} faidx {path_draft_lnk} &&\\
{PATH_BWAMEM} index {path_draft_lnk} &&\\
{PATH_BWAMEM} mem -SP5M -t{n_jobs} {path_draft_lnk} {path_hic_read1} {path_hic_read2} -o {path_output_prefix}.hic_aligned.bam &&\\
{PATH_SAMTOOLS} sort -n -@ {n_jobs} -o {path_output_prefix}.hic_aligned.sorted_name.bam {path_output_prefix}.hic_aligned.bam &&\\
{PATH_JAVA} -jar {PATH_PICARD} MarkDuplicates I={path_output_prefix}.hic_aligned.sorted_name.bam O={path_output_prefix}.hic_aligned.sorted_name.MarkDuplicates.bam M={path_output_prefix}.hic_aligned.MarkDuplicates_metrics.txt
"""

    cmd += f"""
{PATH_YAHS} {path_draft_lnk} {path_output_prefix}.hic_aligned.sorted_name.MarkDuplicates.bam -o {path_output_prefix}.yahs_scaffold"""
    
    for param_name, param_value in additional_params.items():
        if len(param_name) == 1:
            cmd += f" -{param_name} {param_value}"
        else:
            cmd += f" --{param_name} {param_value}"
    
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
        path_draft = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.asm.dip.hap{hap}.p_ctg.fasta"
        path_hic_read1 = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/HiC/Merged/KOREF-S1_PBMC_Hi-C.dedup.trim_F15_T0.merged.2015_DpnII-HinfI_MboI.2019_DpnII-HinfI.R1.fq.gz"
        path_hic_read2 = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/HiC/Merged/KOREF-S1_PBMC_Hi-C.dedup.trim_F15_T0.merged.2015_DpnII-HinfI_MboI.2019_DpnII-HinfI.R2.fq.gz"
        path_output_prefix = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/YaHS.HiC_merged/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.asm.dip.hap{hap}.YaHS.HiC_merged"
        

        modified_option = {
        }

        run(
            path_draft, 
            path_hic_read1, 
            path_hic_read2, 
            path_output_prefix,
            f"{dir_workspace}/YaHS_Scaffolding/tmp",
            f"{dir_workspace}/YaHS_Scaffolding/sh",
            f"{dir_workspace}/YaHS_Scaffolding/log",
            f"YaHS.Hifiasm_Assembly.KOREF_S1.ul_L100kb_QV15.default.hap{hap}.HiC_merged.{curr_time}",
            32,
            [],
            [],
            run_mapping=True,
            **modified_option
        )
#%%