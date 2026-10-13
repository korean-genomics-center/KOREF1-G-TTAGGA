import pandas as pd
import os, time
import subprocess

PYTHON_BIN = "/BiO/Access/yoonsung/anaconda3/bin/python"
JAVA_BIN = "/usr/bin/java"
PATH_CUSTOM_FILTER = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/CallSV_malonge/CallSV/filter.py"
PATH_IRIS_SRC = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/IrisSV/Iris_v1.0.5_fork/Iris/src"
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

def run(path_reference, path_input_vcf, path_input_bam, path_output, rname_first, rname_second, ishifi_first, ishifi_second, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list()):    
    minimap_mode = " --asm5_mixed_racon "
    if ishifi_first:
        minimap_mode += " --hifi_first "
    if ishifi_second:
        minimap_mode += " --hifi_second "

    dir_save_iris_intermediate = os.path.join(os.path.dirname(path_output), "iris_intermediates")

    cmd = f"""{JAVA_BIN} -cp {PATH_IRIS_SRC} Iris {minimap_mode} --keep_long_variants --also_deletions --keep_files minimap_path={PATH_MINIMAP2} samtools_path={PATH_SAMTOOLS} racon_path={PATH_RACON} genome_in={path_reference} vcf_in={path_input_vcf} reads_in={path_input_bam} vcf_out={path_output} out_dir={dir_save_iris_intermediate} threads={n_jobs} rnames_fieldname_fir={rname_first} rnames_fieldname_sec={rname_second}"""
        
    os.makedirs(os.path.dirname(path_output), exist_ok=True)
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
        triobin = dict_hap_to_parental[hap]
            
        path_reference = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Reference_Fasta/Gapclosed_Fasta/CuratedGapclosedASM.KOREF_S1.Trio_Specific.UpperCase.hap{hap}.fasta"
        path_input_bam_first = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Read_Alignment/Align.Winnowmap2/Align_Gapclosed_ASM/Primary_Alignment/Winnowmap_Align.Primary.Ref.KOREF_ASM.hap{hap}.QRY.PacBio_HiFi_ccs_converted.Filter_L1000_QV20.triobin_split.triobin_{triobin}0a.sorted.bam"
        path_input_bam_second = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Read_Alignment/Align.Winnowmap2/Align_Gapclosed_ASM/Primary_Alignment/Winnowmap_Align.Primary.Ref.KOREF_ASM.hap{hap}.QRY.ONT_UL.Filter_L10kb_QV15.triobin_split.triobin_{triobin}0a.sorted.bam"

        path_input_vcf = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Align_Gapclosed_ASM/SV_Call.Inc_Amb/SV_Calling/Merged_SV/SV_Call.KOREF_GapclosedASM_SV_Polishing.Inc_Amb.hap{hap}.hifi_ont.sametrio.filtered.iris.merged.split.homoalt_include.hetalt_only.filt_difftrio_homoref.filt_no_asmsv.vaf_0.5.minsupp_5.call_merged.vcf"
        path_output = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Align_Gapclosed_ASM/SV_Call.Inc_Amb/SV_Calling/Merged_SV/SV_Call.KOREF_GapclosedASM_SV_Polishing.Inc_Amb.hap{hap}.hifi_ont.sametrio.filtered.iris.merged.split.homoalt_include.hetalt_only.filt_difftrio_homoref.filt_no_asmsv.vaf_0.5.minsupp_5.call_merged.iris.vcf"

        jobname_suffix = f"Hap{hap}.ONT_UL.PacBio_HiFi_ccs_converted.triobin_{triobin}0a"

        rname_first = "RNAMES_HIFI"
        rname_second = "RNAMES_ONT"

        assert os.path.exists(path_reference)
        assert os.path.exists(path_input_bam_first)
        assert os.path.exists(path_input_bam_second)
        assert os.path.exists(path_input_vcf)

        curr_time = int(time.time())
        run(
            path_reference,
            path_input_vcf,
            f"{path_input_bam_first},{path_input_bam_second}",
            path_output,
            rname_first,
            rname_second,
            True,
            False,
            f"{dir_workspace}/Polishing/Polishing_SV/IRIS_SV_Filtering/tmp",
            f"{dir_workspace}/Polishing/Polishing_SV/IRIS_SV_Filtering/sh",
            f"{dir_workspace}/Polishing/Polishing_SV/IRIS_SV_Filtering/log",
            f"IRIS_SV_Filtering.{jobname_suffix}.{curr_time}",
            48,
            hostname = [],
            jobs_wait = []
        )