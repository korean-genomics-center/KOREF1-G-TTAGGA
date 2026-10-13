import pandas as pd
import os, time
import subprocess

PATH_WINNNOWMAP = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/winnowmap/winnowmap.v2.03_fork/Winnowmap/bin/winnowmap"
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

def run(dir_annotations, filename_annotjson, filename_biasedname, path_bam, path_save, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list()):    

    cmd = f"""
docker run --ipc=host \\
    -v {dir_annotations}:/annotation \\
    -v {os.path.dirname(path_bam)}:/bam \\
    -v {os.path.dirname(path_save)}:/output \\
    mobinasri/flagger:v1.1.0 bam2cov \\
    --threads {n_jobs} \\
    --baselineAnnotation WHOLE_GENOME --runBiasDetection \\
    --bam /bam/{os.path.basename(path_bam)} \\
    --annotationJson /annotation/{filename_annotjson} \\
    --restrictBiasAnnotations /annotation/{filename_biasedname} \\
    --output /output/{os.path.basename(path_save)}
"""
        
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

    dir_annotations = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Flagger_Evaluation/Annotations/Final_Annotations"
    filename_annotjson = "annotations_path.json"
    filename_biasedname = "biased_annotation_names.txt"

    for datatype in ["ONT_10kb_QV15", "HiFi_1kb_QV20"]:
        path_bam = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Flagger_Evaluation/Winnowmap_Align/Winnowmap_Align.Ref_KOREF_FinalASM_Diploid.Qry_{datatype}.sorted.bam"
        path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Flagger_Evaluation/Bam2Cov/{datatype}/KOREF_FinalASM.Diploid.{datatype}.cov.gz"

        curr_time = int(time.time())
        run(
            dir_annotations, filename_annotjson, filename_biasedname, path_bam, path_save,
            f"{dir_workspace}/Flagger_Evaluation/Run_Bam2Cov/tmp",
            f"{dir_workspace}/Flagger_Evaluation/Run_Bam2Cov/sh",
            f"{dir_workspace}/Flagger_Evaluation/Run_Bam2Cov/log",
            f"Run_Bam2Cov.REF.KOREF_FinalASM.Diploid.Qry.{datatype}.{curr_time}",
            64,
            hostname = [],
            jobs_wait = []
        )