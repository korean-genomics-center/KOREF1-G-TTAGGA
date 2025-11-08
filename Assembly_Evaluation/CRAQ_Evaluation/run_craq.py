import os, time
import subprocess

PATH_CRAQ = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/CRAQ/CRAQ_v1.10/CRAQ/bin/craq"
DIR_MINIMAP2 = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/minimap2/minimap2_2.29_r1283/minimap2-2.29_x64-linux"
DIR_SAMTOOLS = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/Samtools/samtools-1.17/samtools-1.17/bin"


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

def run(path_reference, path_shortread_bam, path_longread_bam, dir_out, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **kwargs):   
    
    cmd = f"""{PATH_CRAQ} -g {path_reference} -ngs {path_shortread_bam} --sms {path_longread_bam} -D {dir_out} -t {n_jobs}"""
    
    for key, val in kwargs.items():
        cmd += f" --{key} {val}"
        
    os.makedirs(os.path.dirname(dir_out), exist_ok=True)
    os.makedirs(dir_tmp, exist_ok=True)
    os.makedirs(dir_shell, exist_ok=True)
    os.makedirs(dir_log, exist_ok=True)
    
    environ = os.environ
    environ["PATH"] = f"{environ['PATH']}:{DIR_MINIMAP2}:{DIR_SAMTOOLS}"
    
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
        ),
        environ = environ
    )



if __name__ == "__main__":
    dict_hap_to_parental = {
        1:'p',
        2:'m'
    }
    
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis"

    for hap in [1, 2]:
        
        path_reference = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Final_Assembly/FinalASM.KOREF_S1.Trio_Specific.hap{hap}.Reordered.fasta"
        
        path_shortread_bam = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/CRAQ_Evaluation/Align.ShortRead.Bwamem2/Bwamem2_Align.Ref.KOREF_FinalASM.hap{hap}.QRY.IlluminaWGS.sorted.bam"
        path_longread_bam = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/CRAQ_Evaluation/Align.LongRead.Winnowmap/Winnowmap_Align.Ref.KOREF_FinalASM.hap{hap}.QRY.ONT_Hifi_Merged.sorted.bam"
        
        dir_out = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/CRAQ_Evaluation/CRAQ_Result/LongRead_Merged_All.ShortRead/Hap{hap}"

        assert os.path.exists(path_reference)
        assert os.path.exists(path_shortread_bam)
        assert os.path.exists(path_longread_bam)

        
        curr_time = int(time.time())
        run(
            path_reference,
            path_shortread_bam,
            path_longread_bam,
            dir_out,
            f"{dir_workspace}/CRAQ_Evaluation/Run_CRAQ/tmp",
            f"{dir_workspace}/CRAQ_Evaluation/Run_CRAQ/sh",
            f"{dir_workspace}/CRAQ_Evaluation/Run_CRAQ/log",
            f"Run_CRAQ.REF.KOREF_FinalASM.hap{hap}.LongRead_Merged_All.ShortRead.{curr_time}",
            64,
            hostname = [],
            jobs_wait = [],
            **{
                "plot" : "T"
            }
        )