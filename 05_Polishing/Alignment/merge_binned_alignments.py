import pandas as pd
import os, time
import subprocess

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

def run(list_alignments, path_save, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list()):    

    cmd = f"{PATH_SAMTOOLS} merge -@{n_jobs} -o {path_save} {' '.join(list_alignments)}"
        
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

    path_align_format = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Read_Alignment/Align.Winnowmap2/Align_PolishedASM.Polish_R8.SV/Primary_Alignment/Winnowmap_Align.Primary.Ref.KOREF_PolishedASM.Polish_R8.SV.hap{hap}.QRY.{tgs_name}.triobin_split.triobin_{triobin}.bam"
    
    curr_time = int(time.time())


    dict_tgs_to_name = {
        "HiFi":"PacBio_HiFi_ccs_converted.Filter_L1000_QV20",
        "ONT":"ONT_UL.Filter_L10kb_QV15"
    }

    dict_hap_to_triobin = {
        1:'p',
        2:'m'
    }
    
    for tgs_type, tgs_name in dict_tgs_to_name.items():
        for hap in [1, 2]:
            
            list_triobinsets = [(dict_hap_to_triobin[hap], '0', 'a')]
            if tgs_type == "HiFi":
                list_triobinsets.append(('p', 'm', '0', 'a'))
            
            
            for triobin_set in list_triobinsets:
                setname = ''.join(triobin_set)

                list_bam = list(map(lambda triobin: path_align_format.format(hap = hap, triobin = triobin, tgs_name = tgs_name), triobin_set))

                path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Read_Alignment/Align.Winnowmap2/Align_PolishedASM.Polish_R8.SV/Primary_Alignment/Winnowmap_Align.Primary.Ref.KOREF_PolishedASM.Polish_R8.SV.hap{hap}.QRY.{tgs_name}.triobin_split.triobin_{setname}.bam"
                
                for path_bam in list_bam:
                    assert os.path.exists(path_bam), path_bam

                curr_time = int(time.time())
                run(
                    list_bam[::-1],
                    path_save,
                    f"{dir_workspace}/Merge_BAM/tmp",
                    f"{dir_workspace}/Merge_BAM/sh",
                    f"{dir_workspace}/Merge_BAM/log",
                    f"Merge_BAM.{tgs_type}.REF.KOREF_PolishedASM.Polish_R8.SV.hap{hap}.triobin_{setname}.{curr_time}",
                    48,
                    hostname = [],
                    jobs_wait = []
                )