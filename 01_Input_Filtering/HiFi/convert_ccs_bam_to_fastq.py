#%%
import subprocess, os, time

PATH_B2F = "/BiO/Access/yoonsung/anaconda3/envs/pbtk/bin/bam2fastq"

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
    script_header += f"#$ -l hostname=Client3|Client4|Client5\n"
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
    
def run_job(path_bam, path_save, jobsuffix, dir_tmp, dir_shell, dir_log, n_jobs = 50):
    cmd = """{PATH_B2F} -o {path_save} --num-threads {n_jobs} {path_bam}""".format(
        PATH_B2F = PATH_B2F,
        path_bam = path_bam,
        path_save = path_save,
        n_jobs = n_jobs
    )
    
    jobname = f"Convert_BAM_to_Fastq_{jobsuffix}.{int(time.time())}"
    
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
            list()
        )
    )


if __name__ == "__main__":
    workspace = "/BiO/Access/yoonsung/Research/KOREF_Assembly/Analysis/Bam2Fastq"
    
    path_bam = "/BiO/Access/yoonsung/Research/KOREF_Assembly/Results/PacBio/Subread_to_CCS/CHM13/m64062_190803_042216.ccs.bam"
    path_save = "/BiO/Access/yoonsung/Research/KOREF_Assembly/Results/PacBio/Subread_to_CCS/CHM13/m64062_190803_042216.ccs.converted"
    
    
    run_job(
        path_bam,
        path_save,
        "CHM13_m64062_190803_042216",
        os.path.join(workspace, "tmp"),
        os.path.join(workspace, "sh"),
        os.path.join(workspace, "log")
    )


if __name__ == "__main2__":
    list_bam_names = [
        "KOREF-S1_Blood_WGS_Sequel-II_20190614_FC00",
        "KOREF-S1_Blood_WGS_Sequel-II_20190709_FC00",
        "KOREF-S1_Blood_WGS_Sequel-II_20190822_FC00",
        "KOREF-S1_Blood_WGS_Sequel-II_20190918_FC00",
        "KOREF-S1_Blood_WGS_Sequel-II_20190923_FC00",
        "KOREF-S1_Blood_WGS_Sequel-II_20191210_FC00"
    ]
    
    workspace = "/BiO/Access/yoonsung/Research/KOREF_Assembly/Analysis/Bam2Fastq"
    
    for name in list_bam_names:
        path_bam = f"/BiO/Access/yoonsung/Research/KOREF_Assembly/Results/PacBio/Subread_to_CCS/{name}/{name}.ccs.bam"
        path_save = f"/BiO/Access/yoonsung/Research/KOREF_Assembly/Results/PacBio/Subread_to_CCS/{name}/{name}.ccs.converted"
    
        run_job(
            path_bam,
            path_save,
            f"KOREFS1_PB_{name}",
            os.path.join(workspace, "tmp"),
            os.path.join(workspace, "sh"),
            os.path.join(workspace, "log")
        )
#%%