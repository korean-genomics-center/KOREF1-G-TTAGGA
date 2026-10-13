#%%
import subprocess, os, time
from copy import deepcopy

DIR_HIFIADAPTERFILT = "/BiO/Access/yoonsung/Research/KOREF_Assembly/Resources/Tool/HiFiAdapterFilt/HiFiAdapterFilt"
DIR_BLAST = "/BiO/Access/yoonsung/Research/KOREF_Assembly/Resources/Tool/HiFiAdapterFilt/BLAST/ncbi-blast-2.16.0+/bin"
DIR_BAMTOOLS = "/BiO/Access/yoonsung/Research/KOREF_Assembly/Resources/Tool/HiFiAdapterFilt/Bamtools/bamtools/build/src"

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
        
    curr_env = deepcopy(os.environ)
    curr_env["PATH"] = f"{curr_env['PATH']}:{DIR_HIFIADAPTERFILT}:{DIR_BLAST}:{DIR_BAMTOOLS}:{DIR_HIFIADAPTERFILT}/DB"
    subprocess.run(f"qsub {path_shell}", shell = True, env=curr_env)
    
def run_job(prefix_input, prefix_save, n_jobs, jobsuffix, dir_tmp, dir_shell, dir_log):
    filename_prefix = os.path.basename(prefix_input)
    dir_file = os.path.dirname(prefix_input)
    os.chdir(dir_file)
    cmd = """cd {dir_file};hifiadapterfilt.sh -p {prefix_input} -o {prefix_save} -t {n_jobs}""".format(
        DIR_HIFIADAPTERFILT = DIR_HIFIADAPTERFILT,
        DIR_BLAST = DIR_BLAST,
        DIR_BAMTOOLS = DIR_BAMTOOLS,
        prefix_input = filename_prefix,
        prefix_save = prefix_save,
        n_jobs = n_jobs,
        dir_file = dir_file
    )
    
    jobname = f"Adapter_Filt_Pacbio_{jobsuffix}.{int(time.time())}"
    
    os.makedirs(prefix_save, exist_ok=True)
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
    workspace = "/BiO/Access/yoonsung/Research/KOREF_Assembly/Analysis/Adapter_Filt_Pacbio"
    
    run_job(
        prefix_input = "/BiO/Access/yoonsung/Research/KOREF_Assembly/Results/PacBio/Subread_to_CCS/CHM13/m64062_190803_042216.ccs.converted",
        prefix_save = f"/BiO/Access/yoonsung/Research/KOREF_Assembly/Results/PacBio/Adapter_Trimmed/CHM13",
        n_jobs = 16,
        jobsuffix = f"CHM13_m64062_190803_042216",
        dir_tmp = os.path.join(workspace, "tmp"),
        dir_shell = os.path.join(workspace, "sh"),
        dir_log = os.path.join(workspace, "log")
    )

if __name__ == "__main2__":
    workspace = "/BiO/Access/yoonsung/Research/KOREF_Assembly/Analysis/Adapter_Filt_Pacbio"
    
    list_bam_names = [
        "KOREF-S1_Blood_WGS_Sequel-II_20190614_FC00",
        "KOREF-S1_Blood_WGS_Sequel-II_20190709_FC00",
        "KOREF-S1_Blood_WGS_Sequel-II_20190822_FC00",
        "KOREF-S1_Blood_WGS_Sequel-II_20190918_FC00",
        "KOREF-S1_Blood_WGS_Sequel-II_20190923_FC00",
        "KOREF-S1_Blood_WGS_Sequel-II_20191210_FC00"
    ]
    
    for name in list_bam_names:
    
        run_job(
            prefix_input = f"/BiO/Access/yoonsung/Research/KOREF_Assembly/Results/PacBio/Subread_to_CCS/{name}/{name}.ccs.converted",
            prefix_save = f"/BiO/Access/yoonsung/Research/KOREF_Assembly/Results/PacBio/Adapter_Trimmed/{name}",
            n_jobs = 16,
            jobsuffix = f"KOREFS1_PB_{name}",
            dir_tmp = os.path.join(workspace, "tmp"),
            dir_shell = os.path.join(workspace, "sh"),
            dir_log = os.path.join(workspace, "log")
        )
# %%
