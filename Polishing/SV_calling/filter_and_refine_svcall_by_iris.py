import pandas as pd
import os, time
import subprocess

PYTHON_BIN = "/BiO/Access/yoonsung/anaconda3/bin/python"
JAVA_BIN = "/usr/bin/java"
PATH_CUSTOM_FILTER = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/CallSV_malonge/CallSV/filter.py"
PATH_IRIS_SRC = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/IrisSV/Iris_v1.0.5/Iris/src"
PATH_MINIMAP2 = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/minimap2/minimap2_2.29_r1283/minimap2-2.29_x64-linux/minimap2"
PATH_SAMTOOLS = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/Samtools/samtools-1.17/samtools-1.17/samtools"
PATH_RACON = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/Racon/Racon_v1.4.3/racon/build/bin/racon"


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

def run(path_reference, path_input_vcf, path_input_bam, path_output_prefix, is_hifi, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list()):    
    minimap_mode = ""
    if is_hifi:
        minimap_mode = " --hifi "

    path_save_filtered_vcf = f"{path_output_prefix}.filtered.vcf"
    path_save_iris_filtered_vcf = f"{path_output_prefix}.filtered.iris.vcf"

    dir_save_iris_intermediate = os.path.join(os.path.dirname(path_output_prefix), "iris_intermediates")

    cmd = f"""{PYTHON_BIN} {PATH_CUSTOM_FILTER} {path_input_vcf} > {path_save_filtered_vcf} && {JAVA_BIN} -cp {PATH_IRIS_SRC} Iris {minimap_mode} --keep_long_variants --also_deletions minimap_path={PATH_MINIMAP2} samtools_path={PATH_SAMTOOLS} racon_path={PATH_RACON} genome_in={path_reference} vcf_in={path_save_filtered_vcf} reads_in={path_input_bam} vcf_out={path_save_iris_filtered_vcf} out_dir={dir_save_iris_intermediate} threads={n_jobs}"""
        
    os.makedirs(os.path.dirname(path_output_prefix), exist_ok=True)
    os.makedirs(dir_save_iris_intermediate, exist_ok=True)
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

    dict_hap_to_parental = {
        1:"p",
        2:"m"
    }

    for hap in [1, 2]:
        for triobins in ['m', 'p']:
            for bam_filter in ["", ".MapQ20.Err0.05"]:
            
                path_reference = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Reference_Fasta/Polished_Fasta/CuratedGapclosedASM.KOREF_S1.Trio_Specific.UpperCase.hap{hap}.fasta"
                path_input_bam = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Read_Alignment/Align.Winnowmap2/Align_Gapclosed_ASM/Primary_Alignment/Winnowmap_Align.Primary.Ref.KOREF_ASM.hap{hap}.QRY.PacBio_HiFi_ccs_converted.Filter_L1000_QV20.triobin_split.triobin_{triobins}.sorted{bam_filter}.bam"
                path_input_vcf = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Align_Gapclosed_ASM/SV_Call/Sniffles2/PacBio_HiFi_ccs_converted.Filter_L1000_QV20/VariantCall.SV_call.Sniffles2.Ref.KOREF_ASM.hap{hap}.Read.PacBio_HiFi_ccs_converted.Filter_L1000_QV20.triobin_{triobins}{bam_filter}.vcf"
                path_output_prefix = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Align_Gapclosed_ASM/SV_Call/Sniffles2/PacBio_HiFi_ccs_converted.Filter_L1000_QV20/Filtered/VariantCall.SV_call.Sniffles2.Ref.KOREF_ASM.hap{hap}.Read.PacBio_HiFi_ccs_converted.Filter_L1000_QV20.triobin_{triobins}{bam_filter}"

                is_hifi = True

                jobname_suffix = f"Hap{hap}.PacBio_HiFi_ccs_converted.Filter_L1000_QV20.triobin_{triobins}{bam_filter}"


                assert os.path.exists(path_reference)
                assert os.path.exists(path_input_bam)
                assert os.path.exists(path_input_vcf)

                curr_time = int(time.time())
                run(
                    path_reference,
                    path_input_vcf,
                    path_input_bam,
                    path_output_prefix,
                    is_hifi,
                    f"{dir_workspace}/Polishing/IRIS_SV_Filtering/tmp",
                    f"{dir_workspace}/Polishing/IRIS_SV_Filtering/sh",
                    f"{dir_workspace}/Polishing/IRIS_SV_Filtering/log",
                    f"IRIS_SV_Filtering.{jobname_suffix}.{curr_time}",
                    48,
                    hostname = [],
                    jobs_wait = []
                )


if __name__ == "__main__":
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis"

    dict_hap_to_parental = {
        1:"p",
        2:"m"
    }

    for hap in [1, 2]:
        for triobins in ['m', 'p']:
            for bam_filter in ["", ".MapQ20.Err0.1"]:
            
                path_reference = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Reference_Fasta/Polished_Fasta/CuratedGapclosedASM.KOREF_S1.Trio_Specific.UpperCase.hap{hap}.fasta"
                path_input_bam = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Read_Alignment/Align.Winnowmap2/Align_Gapclosed_ASM/Primary_Alignment/Winnowmap_Align.Primary.Ref.KOREF_ASM.hap{hap}.QRY.ONT_UL.Filter_L10kb_QV15.triobin_split.triobin_{triobins}.sorted{bam_filter}.bam"
                path_input_vcf = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Align_Gapclosed_ASM/SV_Call/Sniffles2/ONT_UL.Filter_L10kb_QV15/VariantCall.SV_call.Sniffles2.Ref.KOREF_ASM.hap{hap}.Read.ONT_UL.Filter_L10kb_QV15.triobin_{triobins}{bam_filter}.vcf"
                path_output_prefix = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Align_Gapclosed_ASM/SV_Call/Sniffles2/ONT_UL.Filter_L10kb_QV15/Filtered/VariantCall.SV_call.Sniffles2.Ref.KOREF_ASM.hap{hap}.Read.ONT_UL.Filter_L10kb_QV15.triobin_{triobins}{bam_filter}"

                is_hifi = False

                jobname_suffix = f"Hap{hap}.ONT_UL.Filter_L10kb_QV15.triobin_{triobins}{bam_filter}"


                assert os.path.exists(path_reference)
                assert os.path.exists(path_input_bam)
                assert os.path.exists(path_input_vcf)

                curr_time = int(time.time())
                run(
                    path_reference,
                    path_input_vcf,
                    path_input_bam,
                    path_output_prefix,
                    is_hifi,
                    f"{dir_workspace}/Polishing/IRIS_SV_Filtering/tmp",
                    f"{dir_workspace}/Polishing/IRIS_SV_Filtering/sh",
                    f"{dir_workspace}/Polishing/IRIS_SV_Filtering/log",
                    f"IRIS_SV_Filtering.{jobname_suffix}.{curr_time}",
                    48,
                    hostname = [],
                    jobs_wait = []
                )