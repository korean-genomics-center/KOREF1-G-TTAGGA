import os, time
import subprocess

PATH_MINIMAP2 = "/BiO/Research/Project2/Endangered_Jindo_Dog/Resources/Tools/minimap2/minimap2_v2.30_r1287/minimap2/minimap2"


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

def make_symlink(path_src, dir_dst, ch_extension = None):
    filename_dst = os.path.basename(path_src)
    if ch_extension:
        filename_dst = f"{'.'.join(filename_dst.split('.')[:-1])}.{ch_extension}"
    path_dst = os.path.join(dir_dst, filename_dst)
    if os.path.exists(path_dst):
        pass
    else:
        os.symlink(path_src, path_dst)
    return filename_dst

def run(qryfasta, reffasta, path_save, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):
    os.makedirs(os.path.dirname(path_save), exist_ok=True)

    cmd = f"{PATH_MINIMAP2} -x sr -c -N 1000000 --secondary=yes -p 0.5 -o {path_save} -t {n_jobs}"
    
    for param_name, param_value in additional_params.items():
        if len(param_name) == 1:
            cmd += f" -{param_name} {param_value}"
        else:
            cmd += f" --{param_name} {param_value}"
    cmd += f" {reffasta} {qryfasta}"
    
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
    dir_workspace = "/BiO/Research/Project2/Endangered_Jindo_Dog/Analysis"

    curr_time = int(time.time())
    

    for hap in [1, 2]:
        path_ref = f"/BiO/Research/Project2/Endangered_Jindo_Dog/Results/Genome_Assembly/Final_Assembly/FinalASM.J495799_Child.Trio_Specific.hap{hap}.fasta"
        # path_qry = f"/BiO/Research/Project2/Endangered_Jindo_Dog/Results/Repeat_Annotation/Satellite_Repeat/run_srf/Hap{hap}/Satellite_repeat_detection.J495799_Child.Hap{hap}.k151.ci20.srf_output.min_length_100.fa"
        path_qry = "/BiO/Research/Project2/Endangered_Jindo_Dog/Resources/Reference/Repeat_Library/canis_familiaris_clone1_satellite/sequence.fa"
        
        
        # jobname = f"Align.Ref.J495799_Child_Assembly_Hap{hap}.QRY_SRF_result.k151.ci20.min100bp.srcN1000000"

        # path_save = f"/BiO/Research/Project2/Endangered_Jindo_Dog/Results/Repeat_Annotation/Satellite_Repeat/run_srf/Hap{hap}/Repeat_Count/Align.Ref_J495799_Child.Hap{hap}.QRY_SRF_result.k151.ci20.min100bp.srcN1000000.Pcf1000r100100.paf"
        
        jobname = f"Align.Ref.J495799_Child_Assembly_Hap{hap}.QRY_known_satellite.srcN1000000"

        path_save = f"/BiO/Research/Project2/Endangered_Jindo_Dog/Results/Repeat_Annotation/Satellite_Repeat/run_srf/Hap{hap}/Repeat_Count/Align.Ref_J495799_Child.Hap{hap}.QRY_known_satellite.srcN1000000.Pcf1000r100100.paf"

        run(
            path_qry, 
            path_ref, 
            path_save, 
            f"{dir_workspace}/Minimap2_ASM_Align/tmp",
            f"{dir_workspace}/Minimap2_ASM_Align/sh",
            f"{dir_workspace}/Minimap2_ASM_Align/log",
            f"Minimap2.{jobname}.{curr_time}",
            24,
            [],
            []
        )
    