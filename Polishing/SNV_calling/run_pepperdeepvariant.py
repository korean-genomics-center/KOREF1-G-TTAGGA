import os, time
import subprocess


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

def make_symlink(path_src, path_dst, replace = True):
    if os.path.exists(path_dst):
        if replace:
            os.remove(path_dst)
            os.symlink(path_src, path_dst)
    else:
        os.symlink(path_src, path_dst)

def run(dir_input, filename_bam, filename_reference, dir_out, output_name, dir_tmp, dir_shell, dir_log, jobname, n_jobs, pepper_mode, hostname = list(), jobs_wait = list()):    

    assert os.path.exists(os.path.join(dir_input, filename_bam))
    assert os.path.exists(os.path.join(dir_input, f"{filename_bam}.bai"))
    assert os.path.exists(os.path.join(dir_input, filename_reference))
    os.makedirs(dir_out, exist_ok=True)

    cmd = f"""
docker run --ipc=host \\
    -v {dir_input}:/input \\
    -v {dir_out}:/output \\
    -e OMP_NUM_THREADS={n_jobs} \\
    -e TF_NUM_INTRAOP_THREADS={n_jobs} \\
    -e TF_NUM_INTEROP_THREADS=1 \\
    kishwars/pepper_deepvariant:r0.8 \\
    run_pepper_margin_deepvariant call_variant \\
    -b /input/{filename_bam} \\
    -f /input/{filename_reference} \\
    -o /output \\
    -p {output_name} \\
    -t {n_jobs} \\
    --{pepper_mode} --gvcf --phased_output --skip_final_phased_bam
"""
        
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
    
    dict_hap_to_parental = {
        1:"p",
        2:"m"
    }
    
    dict_hap_exchange = {
        'p' : 'm',
        'm' : 'p'
    }
    dict_tgs_to_name = {
        "HiFi":"PacBio_HiFi_ccs_converted.Filter_L1000_QV20",
        "ONT":"ONT_UL.Filter_L10kb_QV15"
    }
    dict_tgs_to_mode = {
        "HiFi":"hifi",
        "ONT":"ont_r10_q20"
    }
    
    curr_time = int(time.time())

    for tgs in ["ONT", "HiFi"]:
        for hap in [1, 2]:
            if tgs == "HiFi":
                list_triobins = ["pm0a", dict_hap_exchange[dict_hap_to_parental[hap]], f"{dict_hap_to_parental[hap]}0a", dict_hap_to_parental[hap]]
            elif tgs == "ONT":
                list_triobins = [dict_hap_exchange[dict_hap_to_parental[hap]], f"{dict_hap_to_parental[hap]}0a", dict_hap_to_parental[hap]]
                # list_triobins = [dict_hap_exchange[dict_hap_to_parental[hap]]]
                
            for triobins in list_triobins:
                
                dir_input = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Read_Alignment/Align.Winnowmap2/Align_PolishedASM.Polish_R3.SV/Primary_Alignment/Docker_Input/{tgs}.Hap{hap}.triobin_{triobins}"

                filename_bam = f"Winnowmap_Align.Primary.Ref.KOREF_PolishedASM.Polish_R3.SV.hap{hap}.QRY.{dict_tgs_to_name[tgs]}.triobin_split.triobin_{triobins}.sorted.bam"
                filename_reference = f"PolishedASM.KOREF_S1.Trio_Specific.UpperCase.Polish_R3.SV.hap{hap}.fasta"
                
                dir_output = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round3_SNV_Call/SNV_Calling/SNV_Calling.{tgs}.Hap{hap}.Triobin_{triobins}"
                
                
                output_name = f"SNV_Call.PepperDeepvariant.KOREF_PolishedASM.Polish_R3.SV.hap{hap}.{tgs}.triobin_{triobins}"

                pepper_mode = dict_tgs_to_mode[tgs]

                jobname_suffix = f"KOREF_PolishedASM.Polish_R3.SV.hap{hap}.{tgs}.triobin_{triobins}"
                
                n_jobs = 48
                
                run(
                    dir_input, 
                    filename_bam, 
                    filename_reference,
                    dir_output,
                    output_name,
                    f"{dir_workspace}/Polishing/Polishing_SNV/Pepper_Variant_Call/tmp",
                    f"{dir_workspace}/Polishing/Polishing_SNV/Pepper_Variant_Call/sh",
                    f"{dir_workspace}/Polishing/Polishing_SNV/Pepper_Variant_Call/log",
                    f"PepperDeepvariant_SNV.{jobname_suffix}.{curr_time}",
                    n_jobs,
                    pepper_mode,
                    hostname = [],
                    jobs_wait = []
                )