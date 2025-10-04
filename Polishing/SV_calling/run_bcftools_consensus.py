import pandas as pd
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

def write_shell_script_and_run(path_shell, shell_context):
    with open(path_shell, 'w') as fw:
        fw.write(shell_context)
    subprocess.run(f"qsub {path_shell}", shell = True)

def run(path_reference, path_vcf, path_save, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list()):   

    cmd = f"""{PATH_BCFTOOLS} annotate -x INFO --force {path_vcf} > {path_vcf}.reheader.vcf && {PATH_BCFTOOLS} sort -Oz {path_vcf}.reheader.vcf > {path_vcf}.reheader.sorted.vcf.gz && {PATH_BCFTOOLS} index {path_vcf}.reheader.sorted.vcf.gz && {PATH_BCFTOOLS} consensus {path_vcf}.reheader.sorted.vcf.gz -f {path_reference} -H A > {path_save}"""
        
    os.makedirs(os.path.dirname(path_save), exist_ok=True)
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

    curr_time = int(time.time())

    for hap in [1, 2]:
        path_reference = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Reference_Fasta/Round3.SV_SNV_Polished_Fasta/PolishedASM.KOREF_S1.Trio_Specific.UpperCase.Polish_R3.SV_SNV.hap{hap}.fasta"
        path_vcf = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round4_SV_Call/SV_Call_Merged/Merged_Call/SV_Call.KOREF_PolishedASM_Round4_SV.hap{hap}.hifi_ont.sametrio.filtered.iris.final_evaluated.indel_only.additional_support_of_unclassified_reads_on_depth_5.manual_curation.filter_overlap.vcf"
        path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Reference_Fasta/Round4.SV_Polished_Fasta/PolishedASM.KOREF_S1.Trio_Specific.UpperCase.Polish_R4.SV.hap{hap}.fasta"

        assert os.path.exists(path_reference)
        assert os.path.exists(path_vcf)

        name = f"Polish_R4.Hap{hap}.SVcall.From_Sniffles2.ONT_HiFi_CallMerged"
        
        curr_time = int(time.time())
        run(
            path_reference,
            path_vcf,
            path_save,
            f"{dir_workspace}/Polishing/Bcftools_Consensus/tmp",
            f"{dir_workspace}/Polishing/Bcftools_Consensus/sh",
            f"{dir_workspace}/Polishing/Bcftools_Consensus/log",
            f"Bcftools_Consensus.{name}.{curr_time}",
            32,
            hostname = [],
            jobs_wait = []
        )