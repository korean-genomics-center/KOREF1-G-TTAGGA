import pandas as pd
import os, time
import subprocess

PATH_SAMTOOLS = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/Samtools/samtools-1.17/samtools-1.17/bin/samtools"


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

def run(path_bam, path_save, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list()):    

    cmd = f"""{PATH_SAMTOOLS} view {path_bam} | cut -f1 |sort | uniq > {path_save}"""
        
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


if __name__ == "__main2__":
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis"
        
    list_region = [
        "chr13:1-1000000",
        "chr13:1-20000000",
        "chrY:20000000",
        "chr15:98753195"
    ]
    
    for hap in ['a', '0']:
        for region in list_region:
            region_name = region.replace(':', '_').replace('-', '_')
        
            path_bam = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/Align.Winnowmap2/Region_filtered/Winnowmap_Align.Ref.CHM13.QRY.PacBio_HiFi_ccs_converted.Filter_L1000_QV20.triobin_split.triobin_{hap}.{region_name}.bam"

            path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/Align.Winnowmap2/Region_filtered/Extract_ReadID_BAM.Ref.CHM13.QRY.PacBio_HiFi_ccs_converted.Filter_L1000_QV20.triobin_split.triobin_{hap}.{region_name}.readid.txt"

            curr_time = int(time.time())
            run(
                path_bam, 
                path_save,
                f"{dir_workspace}/Extract_ReadID_BAM/tmp",
                f"{dir_workspace}/Extract_ReadID_BAM/sh",
                f"{dir_workspace}/Extract_ReadID_BAM/log",
                f"Extract_ReadID_BAM.HiFi.Filter_L1000_QV20.triobin_split.triobin_{hap}.{region_name}.{curr_time}",
                1,
                hostname = [],
                jobs_wait = []
            )


if __name__ == "__main__":
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis"
        
    list_region = [
        "chr13:1-1000000",
        "chr13:1-20000000",
        "chrY:20000000",
        "chr15:98753195"
    ]
    
    for hap in ['a', '0']:
        for region in list_region:
            region_name = region.replace(':', '_').replace('-', '_')
        
            path_bam = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/Align.Winnowmap2/Region_filtered/Winnowmap_Align.Ref.CHM13.QRY.ONT_UL.Filter_L50kb_QV20_Q30Perc80.triobin_split.triobin_{hap}.{region_name}.bam"

            path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/Align.Winnowmap2/Region_filtered/Extract_ReadID_BAM.Ref.CHM13.QRY.ONT_UL.Filter_L50kb_QV20_Q30Perc80.triobin_split.triobin_{hap}.{region_name}.readid.txt"

            curr_time = int(time.time())
            run(
                path_bam, 
                path_save,
                f"{dir_workspace}/Extract_ReadID_BAM/tmp",
                f"{dir_workspace}/Extract_ReadID_BAM/sh",
                f"{dir_workspace}/Extract_ReadID_BAM/log",
                f"Extract_ReadID_BAM.ONT_UL.Filter_L50kb_QV20_Q30Perc80.triobin_split.triobin_{hap}.{region_name}.{curr_time}",
                1,
                hostname = [],
                jobs_wait = []
            )

if __name__ == "__main2__":
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis"
    

    dict_hap_to_parental = {
        1:'p',
        2:'m'
    }
    
    dict_hap_to_list_region = {
        1:[
            "chr13:1-100000",
            # "chr13:1-1000000",
            # "chr13:1-20000000",
            # "chrY:20000000"
        ],
        2:[
            # "chr15:98753195"
        ]
    }
    
    
    for hap in [1, 2]:
        for region in dict_hap_to_list_region[hap]:
            region_name = region.replace(':', '_').replace('-', '_')
        
            path_bam = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/Align.Winnowmap2/Region_filtered/Winnowmap_Align.Ref.CHM13.QRY.PacBio_HiFi_ccs_converted.Filter_L1000_QV20.triobin_split.triobin_{dict_hap_to_parental[hap]}.{region_name}.bam"

            path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/Align.Winnowmap2/Region_filtered/Extract_ReadID_BAM.Ref.CHM13.QRY.PacBio_HiFi_ccs_converted.Filter_L1000_QV20.triobin_split.triobin_{dict_hap_to_parental[hap]}.{region_name}.readid.txt"

            curr_time = int(time.time())
            run(
                path_bam, 
                path_save,
                f"{dir_workspace}/Extract_ReadID_BAM/tmp",
                f"{dir_workspace}/Extract_ReadID_BAM/sh",
                f"{dir_workspace}/Extract_ReadID_BAM/log",
                f"Extract_ReadID_BAM.HiFi.Filter_L1000_QV20.triobin_split.triobin_{dict_hap_to_parental[hap]}.{region_name}.{curr_time}",
                8,
                hostname = [],
                jobs_wait = []
            )



if __name__ == "__main2__":
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis"
    

    dict_hap_to_parental = {
        1:'p',
        2:'m'
    }
    
    dict_hap_to_list_region = {
        1:[
            "chr13:1-100000",
            # "chr13:1-1000000",
            # "chr13:1-20000000",
            # "chrY:20000000"
        ],
        2:[
            # "chr15:98753195"
        ]
    }
    
    
    for hap in [1, 2]:
        for region in dict_hap_to_list_region[hap]:
            region_name = region.replace(':', '_').replace('-', '_')
        
            path_bam = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/Align.Winnowmap2/Region_filtered/Winnowmap_Align.Ref.CHM13.QRY.ONT_UL.Filter_L50kb_QV20_Q30Perc80.triobin_split.triobin_{dict_hap_to_parental[hap]}.{region_name}.bam"

            path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/Align.Winnowmap2/Region_filtered/Extract_ReadID_BAM.Ref.CHM13.QRY.ONT_UL.Filter_L50kb_QV20_Q30Perc80.triobin_split.triobin_{dict_hap_to_parental[hap]}.{region_name}.readid.txt"

            curr_time = int(time.time())
            run(
                path_bam, 
                path_save,
                f"{dir_workspace}/Extract_ReadID_BAM/tmp",
                f"{dir_workspace}/Extract_ReadID_BAM/sh",
                f"{dir_workspace}/Extract_ReadID_BAM/log",
                f"Extract_ReadID_BAM.ONT_UL.Filter_L50kb_QV20_Q30Perc80.triobin_split.triobin_{dict_hap_to_parental[hap]}.{region_name}.{curr_time}",
                8,
                hostname = [],
                jobs_wait = []
            )
#%%