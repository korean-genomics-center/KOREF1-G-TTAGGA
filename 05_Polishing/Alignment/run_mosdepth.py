import pandas as pd
import os, time
import subprocess

PATH_MOSDEPTH = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/mosdepth/mosdepth_v0.3.11/mosdepth"


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

def run(path_bam, path_reference, path_output_prefix, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list()):    

    cmd = f"""{PATH_MOSDEPTH} -t {n_jobs} --fasta {path_reference} {path_output_prefix} {path_bam}"""
        
    os.makedirs(os.path.dirname(path_output_prefix), exist_ok=True)
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

    dict_tgs_to_name = {
        "ONT": "ONT_UL.Filter_L10kb_QV15",
        "HiFi": "PacBio_HiFi_ccs_converted.Filter_L1000_QV20"
    }

    dict_hap_to_parental = {
        1:'p',
        2:'m'
    }

    for tgs, name in dict_tgs_to_name.items():
        for hap in [1, 2]:
            triobin = dict_hap_to_parental[hap]
                
            path_fasta_ref = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Reference_Fasta/Manual_Polishing/After_Round4/Round7.Manual_Polished_Fasta/PolishedASM.KOREF_S1.Trio_Specific.UpperCase.Polish_R7.Man.hap{hap}.fasta"

            path_bam = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Read_Alignment/Align.Winnowmap2/Align_PolishedASM.Polish_R7.Man/Primary_Alignment/Sorted_Bam/Winnowmap_Align.Primary.Ref.KOREF_PolishedASM.Polish_R7.Man.hap{hap}.QRY.{name}.triobin_split.triobin_{triobin}.sorted.bam"

            path_output_prefix = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round8_SV_Call/SV_Call_without_Unclassified/Alignment_Depth/Mosdepth.Primary.Ref.KOREF_PolishedASM.Polish_R7.Man.hap{hap}.QRY.{name}.triobin_split.triobin_{triobin}"

            assert os.path.exists(path_fasta_ref)
            assert os.path.exists(path_bam)

            curr_time = int(time.time())
            run(
                path_bam,
                path_fasta_ref,
                path_output_prefix,
                f"{dir_workspace}/Polishing/Mosdepth_Depth/tmp",
                f"{dir_workspace}/Polishing/Mosdepth_Depth/sh",
                f"{dir_workspace}/Polishing/Mosdepth_Depth/log",
                f"Mosdepth_Depth.{tgs}.REF.KOREF_PolishedASM.Polish_R7.Man.hap{hap}.triobin_{triobin}.{curr_time}",
                12,
                hostname = [],
                jobs_wait = []
            )