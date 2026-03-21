#%%
import subprocess, os, time

PATH_SALMON = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/Salmon/Salmon_v1.10.0/salmon-latest_linux_x86_64/bin/salmon"

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

def write_shell_script_and_run(path_shell, shell_context, env = None):
    with open(path_shell, 'w') as fw:
        fw.write(shell_context)
    
    subprocess.run(f"qsub {path_shell}", shell = True, env = env)
    # print(path_shell)

def run(list_rnaseq_read1, list_rnaseq_read2, path_salmon_index, dir_output, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):
    # dir_output = os.path.dirname(path_output_prefix)

    cmd = f"""
{PATH_SALMON} quant -i {path_salmon_index} -l A -1 {' '.join(list_rnaseq_read1)} -2 {' '.join(list_rnaseq_read2)} --validateMappings -o {dir_output} --threads {n_jobs}
"""
    
    os.makedirs(dir_output, exist_ok=True)
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

    curr_time = int(time.time())

    dict_rnaseq_data = {
        "Blood_NovaSeq_20190403":{
            "read1":[
                "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WTS/RemoveAdapter/Blood/KOREF-S1_Blood_WTS_NovaSeq_20190403_L01_R1.trimPolyG.Trim_F15_T5.fq.gz"
            ],
            "read2":[
                "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WTS/RemoveAdapter/Blood/KOREF-S1_Blood_WTS_NovaSeq_20190403_L01_R2.trimPolyG.Trim_F15_T5.fq.gz"
            ]
        },
        "CellLine_HiSeq_20180417":{
            "read1":[
                "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WTS/RemoveAdapter/CellLine/KOREF-S1_CellLine_WTS_HiSeq_20180417_L01_R1.trimPolyG.Trim_F15_T5.fq.gz"
            ],
            "read2":[
                "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WTS/RemoveAdapter/CellLine/KOREF-S1_CellLine_WTS_HiSeq_20180417_L01_R2.trimPolyG.Trim_F15_T5.fq.gz"
            ]
        }
    }

    for hap in [1, 2]:
        path_salmon_index = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/CDS_Fasta/Hap{hap}.Blood_Cellline.CGP_PB/Salmon/Index_Salmon/KOREF1_ASM_Hap{hap}.whole_gff.salmon_index"

        for rna_type in dict_rnaseq_data.keys():
            dir_output = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/CDS_Fasta/Hap{hap}.Blood_Cellline.CGP_PB/Salmon/Quant_Salmon/RNASeq_{rna_type}"
            
            run(
                dict_rnaseq_data[rna_type]["read1"], dict_rnaseq_data[rna_type]["read2"], path_salmon_index, dir_output,
                f"{dir_workspace}/CAT_Annotation/Run_Salmon_Quant/tmp",
                f"{dir_workspace}/CAT_Annotation/Run_Salmon_Quant/sh",
                f"{dir_workspace}/CAT_Annotation/Run_Salmon_Quant/log",
                f"Run_Salmon_Quant.hap{hap}.{rna_type}.{curr_time}",
                64,
                [],
                []
            )