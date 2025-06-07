import os, time
import subprocess

PATH_UNIMAP = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/unimap/unimap/unimap"


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

    cmd = f"{PATH_UNIMAP} -x asm5 -o {path_save} -t {n_jobs}"
    
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

def run_between_two_directories(dir_ref, dir_qry, path_output_prefix, dir_tmp, dir_shell, dir_log, jobname_prefix, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):
    import os
    curr_time = int(time.time())

    list_files_ref = os.listdir(dir_ref)
    list_files_qry = os.listdir(dir_qry)

    for file_ref in list_files_ref:
        for file_qry in list_files_qry:
            prefix_ref = file_ref.replace(".fasta.gz", '').replace(".fa.gz", '').replace(".fasta", '').replace(".fa", '')
            prefix_qry = file_qry.replace(".fasta.gz", '').replace(".fa.gz", '').replace(".fasta", '').replace(".fa", '')
            
            path_output = f"{path_output_prefix}.Ref_{prefix_ref}.Qry_{prefix_qry}.paf"

            run(
                os.path.join(dir_qry, file_qry), 
                os.path.join(dir_ref, file_ref), 
                path_output, 
                dir_tmp,
                dir_shell,
                dir_log,
                f"{jobname_prefix}.Ref_{prefix_ref}.Qry_{prefix_qry}.{curr_time}",
                n_jobs,
                [],
                [],
                **additional_params
            )

if __name__ == "__main__":
    hap = 1

    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis"
    
    dir_ref_seqs = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling.RagTag.default_option.CHM13/DraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.telo_single.hap{hap}.RagTag.default_option.minimap2_latest.CHM13/OneToOneMapping/Raw_Fasta/Reference"
    dir_qry_seqs = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling.RagTag.default_option.CHM13/DraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.telo_single.hap{hap}.RagTag.default_option.minimap2_latest.CHM13/OneToOneMapping/Raw_Fasta/Scaffolds"
    
    path_output_prefix = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling.RagTag.default_option.CHM13/DraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.telo_single.hap{hap}.RagTag.default_option.minimap2_latest.CHM13/OneToOneMapping/Alignment.Unimap/Align.Between_Unmatched.ASM5"

    curr_time = int(time.time())

    modified_option = {
    }

    run_between_two_directories(
        dir_ref_seqs, 
        dir_qry_seqs, 
        path_output_prefix, 
        f"{dir_workspace}/Unimap_ASM_Align/tmp",
        f"{dir_workspace}/Unimap_ASM_Align/sh",
        f"{dir_workspace}/Unimap_ASM_Align/log",
        "Align.Between_Unmatched", 
        16, 
        hostname = list(), 
        jobs_wait = list(), 
        **modified_option)


#%%