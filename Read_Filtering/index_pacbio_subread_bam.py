#%%
import subprocess, os, time

PATH_PBINDEX = "/BiO/Access/yoonsung/anaconda3/envs/ccs/bin/pbindex"

def get_context_for_qsub_submission(cmd, jobname, stderr_path, stdout_path, nthreads, hold_jids = list()):
    #Header
    script_header = f"#!bin/sh\n"
    script_header += f"#$ -N {jobname}\n"
    script_header += f"#$ -S /bin/bash\n"
    script_header += f"#$ -V\n"
    if len(hold_jids) > 0:
        script_header += f"#$ -hold_jid {','.join(hold_jids)}\n"
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
    
def run_job(path_bam, jobsuffix, dir_tmp, dir_shell, dir_log):
    cmd = """{PATH_PBINDEX} {path_bam}""".format(
        PATH_PBINDEX = PATH_PBINDEX,
        path_bam = path_bam
    )
    
    jobname = f"Index_Pacbio_Bam_{jobsuffix}.{int(time.time())}"
    
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
            5,
            list()
        )
    )


# if __name__ == "__main__":
#     workspace = "/BiO/Access/yoonsung/Research/Test_m6A_Methylation/Analysis/Subread_to_CCS"
    
#     path_bam = "/BiO/Access/yoonsung/Research/Test_m6A_Methylation/Resources/RawData/KOREF-S1_Blood_WTS_Sequel-II_20230519_FC00.subreads.bam"
#     path_save = "/BiO/Access/yoonsung/Research/Test_m6A_Methylation/Results/Subread_to_CCS/KOREF-S1_Blood_WTS_Sequel-II_20230519_FC00.ccs.bam"
    
    
#     run_job(
#         path_bam,
#         path_save,
#         "KOREFS1_SequelII_20230519_FC00",
#         os.path.join(workspace, "tmp"),
#         os.path.join(workspace, "sh"),
#         os.path.join(workspace, "log")
#     )


if __name__ == "__main__":
    list_subread_bam = [
        "/BiO/Store/KOREF/KOREF-S1/Blood/WGS/Sequel-II/20190614/KOREF-S1_Blood_WGS_Sequel-II_20190614_FC00.subreads.bam",
        "/BiO/Store/KOREF/KOREF-S1/Blood/WGS/Sequel-II/20190709/KOREF-S1_Blood_WGS_Sequel-II_20190709_FC00.subreads.bam",
        "/BiO/Store/KOREF/KOREF-S1/Blood/WGS/Sequel-II/20190822/KOREF-S1_Blood_WGS_Sequel-II_20190822_FC00.subreads.bam",
        "/BiO/Store/KOREF/KOREF-S1/Blood/WGS/Sequel-II/20190918/KOREF-S1_Blood_WGS_Sequel-II_20190918_FC00.subreads.bam",
        "/BiO/Store/KOREF/KOREF-S1/Blood/WGS/Sequel-II/20190923/KOREF-S1_Blood_WGS_Sequel-II_20190923_FC00.subreads.bam",
        "/BiO/Store/KOREF/KOREF-S1/Blood/WGS/Sequel-II/20191210/KOREF-S1_Blood_WGS_Sequel-II_20191210_FC00.subreads.bam"
    ]
    
    workspace = "/BiO/Access/yoonsung/Research/Test_m6A_Methylation/Analysis/Subread_to_CCS_WGS"
    
    for path_bam in list_subread_bam:
        name = os.path.basename(path_bam).split('.')[0]
    
        run_job(
            path_bam,
            f"KOREFS1_PB_{name}",
            os.path.join(workspace, "tmp"),
            os.path.join(workspace, "sh"),
            os.path.join(workspace, "log")
        )
# %%
