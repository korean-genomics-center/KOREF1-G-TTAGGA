import pandas as pd
import os, time
import subprocess

PATH_JASMINE = "/BiO/Access/yoonsung/anaconda3/envs/jasmine/bin/jasmine"
PATH_SAMTOOLS = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/Samtools/samtools-1.17/samtools-1.17/samtools"


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

def run(list_vcf_files, path_genome, path_output_prefix, dir_tmp, dir_shell, dir_log, jobname, n_jobs, min_supp = 3, hostname = list(), jobs_wait = list()):    
    os.makedirs(os.path.dirname(path_output_prefix), exist_ok=True)

    path_vcf_filelist = f"{path_output_prefix}.input_vcf_list.txt"
    path_save_merged_vcf = f"{path_output_prefix}.merged.vcf"
    dir_save_intermediates = os.path.join(os.path.dirname(path_output_prefix), "jasmine_intermediates")

    with open(path_vcf_filelist, 'w') as fw:
        for path_vcf in list_vcf_files:
            fw.write(f"{path_vcf}\n")

    # max_dist=500 
    cmd = f"""{PATH_JASMINE} file_list={path_vcf_filelist} out_file={path_save_merged_vcf} genome_file={path_genome} max_dist=500 min_seq_id=0.3 spec_reads={min_supp} --dup_to_ins threads={n_jobs} --output_genotypes samtools_path={PATH_SAMTOOLS} out_dir={dir_save_intermediates}"""
        
    os.makedirs(dir_save_intermediates, exist_ok=True)
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
        1:'p',
        2:'m'
    }
    dict_parental_exchange = {
        'p':'m',
        'm':'p'
    }

    for hap in [1, 2]:
        list_vcf_files = [
            f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Align_Gapclosed_ASM/SV_Call/Sniffles2/ONT_UL.Filter_L10kb_QV15/Filtered/VariantCall.SV_call.Sniffles2.Ref.KOREF_ASM.hap{hap}.Read.ONT_UL.Filter_L10kb_QV15.triobin_{dict_hap_to_parental[hap]}.filtered.iris.vcf",
            f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Align_Gapclosed_ASM/SV_Call/Sniffles2/PacBio_HiFi_ccs_converted.Filter_L1000_QV20/Filtered/VariantCall.SV_call.Sniffles2.Ref.KOREF_ASM.hap{hap}.Read.PacBio_HiFi_ccs_converted.Filter_L1000_QV20.triobin_{dict_hap_to_parental[hap]}.filtered.iris.vcf"
        ]

        path_output_prefix = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Align_Gapclosed_ASM/SV_Call/Sniffles2/ONT_HiFi_Merged/VariantCall.SV_call.Sniffles2.Ref.KOREF_ASM.hap{hap}.Read.ONT_UL.Filter_L10kb_QV15.PacBio_HiFi_ccs_converted.Filter_L1000_QV20.triobin_{dict_hap_to_parental[hap]}.filtered.iris"

        jobname_suffix = f"merge_ONT_UL_PacBio_HiFi_ccs_converted_triobin_{dict_hap_to_parental[hap]}"

        for path_vcf in list_vcf_files:
            assert os.path.exists(path_vcf), path_vcf

        curr_time = int(time.time())
        run(
            list_vcf_files, 
            path_output_prefix,
            f"{dir_workspace}/Polishing/Jasmine_Merge_SVcall/tmp",
            f"{dir_workspace}/Polishing/Jasmine_Merge_SVcall/sh",
            f"{dir_workspace}/Polishing/Jasmine_Merge_SVcall/log",
            f"Jasmine_Merge_SVcall.Hap{hap}.{jobname_suffix}.{curr_time}",
            8,
            hostname = [],
            jobs_wait = []
        )