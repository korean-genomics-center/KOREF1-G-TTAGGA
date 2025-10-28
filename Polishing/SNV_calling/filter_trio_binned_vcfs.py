import os, time
import subprocess

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

def run(path_input, path_reference, path_output_prefix, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list()):    
    cmd = f"""
{PATH_BCFTOOLS} norm -m - -f {path_reference} -Oz {path_input} > {path_output_prefix}.Multiallelic_Norm.vcf.gz;
{PATH_BCFTOOLS} view -i 'FORMAT/VAF > 0.5' --threads {n_jobs} -Oz {path_output_prefix}.Multiallelic_Norm.vcf.gz > {path_output_prefix}.Multiallelic_Norm.VAF05.vcf.gz;
{PATH_BCFTOOLS} view --exclude 'GT="RR"' --threads {n_jobs} -Oz {path_output_prefix}.Multiallelic_Norm.VAF05.vcf.gz > {path_output_prefix}.Multiallelic_Norm.VAF05.NonHomRef.vcf.gz

"""
        
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

    for tgs in ["HiFi", "ONT"]:
        for hap in [1, 2]:
            for triobins in [f"{dict_hap_to_parental[hap]}0a", dict_hap_to_parental[hap], dict_hap_exchange[dict_hap_to_parental[hap]]]:
                
                path_input = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round8_SNV_Call/SNV_Filtering/Pass_filter/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.Polish_R8.SV.hap{hap}.{tgs}.triobin_{triobins}.PASS.NonHomRef.vcf.gz"
                path_reference = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Reference_Fasta/Round8.SV_Polished_Fasta/PolishedASM.KOREF_S1.Trio_Specific.UpperCase.Polish_R8.SV.hap{hap}.fasta"

                path_output_prefix = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round8_SNV_Call/SNV_Filtering/Triobinned_for_Evaluation/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.Polish_R8.SV.hap{hap}.{tgs}.triobin_{triobins}.PASS.NonHomRef"

                assert os.path.exists(path_input)
                assert os.path.exists(path_reference)

                jobname_suffix = f"KOREF_PolishedASM.Polish_R8.SV.hap{hap}.{tgs}.triobin_{triobins}.Multiallelic_Norm.VAF05"

                curr_time = int(time.time())
                run(
                    path_input,
                    path_reference,
                    path_output_prefix,
                    f"{dir_workspace}/Polishing/Polishing_SNV/Filter_VCF/tmp",
                    f"{dir_workspace}/Polishing/Polishing_SNV/Filter_VCF/sh",
                    f"{dir_workspace}/Polishing/Polishing_SNV/Filter_VCF/log",
                    f"Filter_VCF.{jobname_suffix}.{curr_time}",
                    8,
                    hostname = [],
                    jobs_wait = []
                )