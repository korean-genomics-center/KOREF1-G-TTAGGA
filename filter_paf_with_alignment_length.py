import os, time
import subprocess

PATH_RAGTAG = "/BiO/Access/yoonsung/anaconda3/envs/RagTag/bin/ragtag.py"
PATH_MUMMERPLOT = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/Mummer/Mummer_release4.0.1/mummer/bin/bin/mummerplot"


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

def run(path_paf, path_paf_filtered, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list()):

    cmd = f"awk '$10 >= 10000 && $12 >= 60' {path_paf} > {path_paf_filtered}"
    
    os.makedirs(os.path.dirname(path_paf_filtered), exist_ok=True)
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

    list_path_qry_format = [
        (
            "KOREF_S1.ul_L100kb_QV15.Trio_Specific.hap{hap}.telo_both",
            "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/Telomere_Detection/Fasta_Split_by_Telomere/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg.telo_both.fasta",
            "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/Telomere_Detection/Fasta_Split_by_Telomere/Minimap2_Align/Align.KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg.telo_both"
        ),
        (
            "KOREF_S1.ul_L100kb_QV15.default_option.hap{hap}.telo_both",
            "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/Telomere_Detection/Fasta_Split_by_Telomere/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.asm.dip.hap{hap}.p_ctg.telo_both.fasta",
            "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/Telomere_Detection/Fasta_Split_by_Telomere/Minimap2_Align/Align.KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.asm.dip.hap{hap}.p_ctg.telo_both"
        )
    ]
    
    dict_refname_to_path = {
        "CHM13" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/CHM13/chm13v2.0.fa",
        "CN1_combine" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/CN1/CN1_combine.v1.0.1.fasta"
    }
    
    curr_time = int(time.time())

    modified_option = {
    }

    for refname, path_ref in dict_refname_to_path.items():
        for hap in ["1", "2"]:
            for jobname_format, path_qry_format, path_save_prefix_format in list_path_qry_format:
                path_save_prefix = path_save_prefix_format.format(hap = hap)
                jobname = f"{jobname_format.format(hap = hap)}.Ref_{refname}"

                path_paf = f"{path_save_prefix}.Ref_{refname}.paf"
                path_save_filtered_paf = f"{path_save_prefix}.Ref_{refname}.Filt_Match10kb_MQ60.paf"

                run(
                    path_paf,
                    path_save_filtered_paf,
                    f"{dir_workspace}/PAF_Filtering/tmp",
                    f"{dir_workspace}/PAF_Filtering/sh",
                    f"{dir_workspace}/PAF_Filtering/log",
                    f"Filter_PAF.{jobname}.{curr_time}",
                    1,
                    [],
                    [],
                    **modified_option
                )
#%%