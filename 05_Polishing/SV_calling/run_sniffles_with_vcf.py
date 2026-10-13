import os, time
import subprocess

PATH_SNIFFLES = "/BiO/Access/yoonsung/anaconda3/envs/sniffles/bin/sniffles"


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

def run(path_input_bam, path_reference, path_genotype, path_save_vcf, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list()):    

    cmd = f"""{PATH_SNIFFLES} -i {path_input_bam} -v {path_save_vcf} --reference {path_reference} -t {n_jobs} --tmp-dir {dir_tmp} --minsupport 3 --mapq 20 --output-rnames --allow-overwrite --genotype-vcf {path_genotype}"""
        
    os.makedirs(os.path.dirname(path_save_vcf), exist_ok=True)
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
    
    dict_hap_exchange = {
        'p' : 'm',
        'm' : 'p'
    }

    dict_bam_format = {
        "PacBio_HiFi_ccs_converted" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Read_Alignment/Align.Winnowmap2/Align_Gapclosed_ASM/Primary_Alignment/Winnowmap_Align.Primary.Ref.KOREF_ASM.hap{hap}.QRY.PacBio_HiFi_ccs_converted.Filter_L1000_QV20.triobin_split.triobin_{triobins}.sorted.MapQ20.Err0.05.bam",
        "ONT_UL" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Read_Alignment/Align.Winnowmap2/Align_Gapclosed_ASM/Primary_Alignment/Winnowmap_Align.Primary.Ref.KOREF_ASM.hap{hap}.QRY.ONT_UL.Filter_L10kb_QV15.triobin_split.triobin_{triobins}.sorted.MapQ20.Err0.1.bam"
    }


    for hap in [1, 2]:
        triobin_opp = dict_hap_exchange[dict_hap_to_parental[hap]]
        for tgs_type in dict_bam_format.keys():

            path_input_bam = dict_bam_format[tgs_type].format(triobins = triobin_opp, hap = hap)

            path_reference = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Reference_Fasta/Polished_Fasta/CuratedGapclosedASM.KOREF_S1.Trio_Specific.UpperCase.hap{hap}.fasta"

            path_genotype = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Align_Gapclosed_ASM/SV_Call/Sniffles2/ONT_HiFi_Merged/VariantCall.SV_call.Sniffles2.Ref.KOREF_ASM.hap{hap}.Read.ONT_UL.Filter_L10kb_QV15.PacBio_HiFi_ccs_converted.Filter_L1000_QV20.triobin_{dict_hap_to_parental[hap]}.MapQ20.Err0.1.filtered.iris.merged.split.hetalt_only.vcf"

            path_save_vcf = path_genotype.replace(".vcf", f".opp_trio.{tgs_type}.vcf")

            jobname_suffix = f"Hap{hap}.{tgs_type}.triobin_{triobin_opp}.Genotype_ONT_UL_Pacbio_Merged_triobin_{dict_hap_to_parental[hap]}"


            assert os.path.exists(path_reference)
            assert os.path.exists(path_input_bam)
            assert os.path.exists(path_genotype)

            curr_time = int(time.time())
            run(
                path_input_bam, 
                path_reference, 
                path_genotype,
                path_save_vcf,
                f"{dir_workspace}/Polishing/Sniffles_SV/tmp",
                f"{dir_workspace}/Polishing/Sniffles_SV/sh",
                f"{dir_workspace}/Polishing/Sniffles_SV/log",
                f"Sniffles_SV.{jobname_suffix}.{curr_time}",
                24,
                hostname = [],
                jobs_wait = []
            )