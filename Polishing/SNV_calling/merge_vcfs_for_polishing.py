import os, time
import subprocess

PATH_BGZIP = "/BiO/Access/yoonsung/anaconda3/envs/bcftools/bin/bgzip"
PATH_BCFTOOLS = "/BiO/Access/yoonsung/anaconda3/envs/bcftools/bin/bcftools"


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

def write_shell_script_and_run(path_shell, shell_context, environ = None):
    with open(path_shell, 'w') as fw:
        fw.write(shell_context)
    subprocess.run(f"qsub {path_shell}", shell = True, env = environ)

def make_symlink(path_src, path_dst, replace = True):
    if os.path.exists(path_dst):
        if replace:
            os.remove(path_dst)
            os.symlink(path_src, path_dst)
    else:
        os.symlink(path_src, path_dst)

def run(list_input, path_output_prefix, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list()):    
    cmd = ""
    
    list_gzipped_input = list()
    for path_input in list_input:
        if not path_input.endswith(".gz"):
            cmd += f"{PATH_BGZIP} -c {path_input} > {path_input}.gz\n"
            path_input = f"{path_input}.gz"
        cmd += f"{PATH_BCFTOOLS} index -f --threads {n_jobs} {path_input}\n"
        list_gzipped_input.append(path_input)
    cmd += f"""
{PATH_BCFTOOLS} concat -a --threads {n_jobs} -Oz {' '.join(list_gzipped_input)} > {path_output_prefix}.vcf.gz;
{PATH_BCFTOOLS} sort --write-index -Oz -o {path_output_prefix}.sorted.vcf.gz {path_output_prefix}.vcf.gz"""
        
    os.makedirs(dir_tmp, exist_ok=True)
    os.makedirs(dir_shell, exist_ok=True)
    os.makedirs(dir_log, exist_ok=True)
    os.makedirs(os.path.dirname(path_output_prefix), exist_ok=True)
    
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
    
    curr_time = int(time.time())

    for hap in [1, 2]:
            
        list_vcfs = [
            f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Align_Polished_ASM.From_R1.SV/SNV_Call/SNV_Calling.HiFi.Hap{hap}.Triobin_pm0a/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.R1_Sv.hap{hap}.HiFi.triobin_pm0a.PASS.NonHomRef.biallelic.GQ30.VAF04.Norm.HomAlt.vcf.gz",
            f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Align_Polished_ASM.From_R1.SV/SNV_Call/SNV_Calling.HiFi.Hap{hap}.Triobin_pm0a/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.R1_Sv.hap{hap}.HiFi.triobin_pm0a.PASS.NonHomRef.biallelic.GQ30.VAF04.Norm.Het.ForPolishing.vcf",
            f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Align_Polished_ASM.From_R1.SV/SNV_Call/SNV_Calling.HiFi.Hap{hap}.Triobin_pm0a/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.R1_Sv.hap{hap}.HiFi.triobin_pm0a.PASS.NonHomRef.multiallelic.GQ25.MinVAF03.SumVAF07.ID_Annot.Norm.ForPolishing.vcf"
        ]
        
        path_output_prefix = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Align_Polished_ASM.From_R1.SV/Polishing_Target/Merged_VCF/HiFi.Hap{hap}.Triobin_pm0a/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.R1_Sv.hap{hap}.HiFi.triobin_pm0a.PASS.NonHomRef.both_biallelic_multiallelic.biallelic_VAF04.ForPolishing"

        for path in list_vcfs:
            assert os.path.exists(path)

        jobname_suffix = f"KOREF_PolishedASM.R1_Sv.hap{hap}"

        run(
            list_vcfs,
            path_output_prefix,
            f"{dir_workspace}/Polishing/Polishing_SNV/Merge_VCF/tmp",
            f"{dir_workspace}/Polishing/Polishing_SNV/Merge_VCF/sh",
            f"{dir_workspace}/Polishing/Polishing_SNV/Merge_VCF/log",
            f"Merge_VCF.{jobname_suffix}.{curr_time}",
            8,
            hostname = [],
            jobs_wait = []
        )