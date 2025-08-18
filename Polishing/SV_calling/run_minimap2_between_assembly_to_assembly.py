import os, time
import subprocess

PATH_MINIMAP2 = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/minimap2/minimap2_2.29_r1283/minimap2-2.29_x64-linux/minimap2"
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

    cmd = f"{PATH_MINIMAP2} -cx asm5 --cs -t {n_jobs}"
    
    for param_name, param_value in additional_params.items():
        if len(param_name) == 1:
            cmd += f" -{param_name} {param_value}"
        else:
            cmd += f" --{param_name} {param_value}"
    cmd += f" {reffasta} {qryfasta} > {path_save}"
    # cmd += f" {reffasta} {qryfasta} | {PATH_SAMTOOLS} view -bS - > {path_save}"
    
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

    curr_time = int(time.time())

    modified_option = {
    }

    dict_ref = {
        "hap1" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Reference_Fasta/Polished_Fasta/CuratedGapclosedASM.KOREF_S1.Trio_Specific.UpperCase.hap1.fasta",
        "hap2" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Reference_Fasta/Polished_Fasta/CuratedGapclosedASM.KOREF_S1.Trio_Specific.UpperCase.hap2.fasta"
    }

    path_output_format = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Align_Gapclosed_ASM/Between_Haplotype/Before_SV_Polishing/Ref_{refname}.Qry_{qryname}/Align_Minimap2.Ref_{refname}.Qry_{qryname}.sorted.paf"

    for refname in dict_ref.keys():
        for qryname in dict_ref.keys():
            if refname == qryname:
                continue

            jobname = f"Ref_{refname}.Qry_{qryname}"
            
            run(
                dict_ref[qryname], 
                dict_ref[refname], 
                path_output_format.format(refname = refname, qryname = qryname), 
                f"{dir_workspace}/Polishing/Minimap2_ASM_Align/tmp",
                f"{dir_workspace}/Polishing/Minimap2_ASM_Align/sh",
                f"{dir_workspace}/Polishing/Minimap2_ASM_Align/log",
                f"Minimap2.{jobname}.{curr_time}",
                48,
                [],
                [],
                **modified_option
            )
            
if __name__ == "__main__":
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis"

    curr_time = int(time.time())

    modified_option = {
    }

    dict_ref_bef = {
        "hap1" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon/Final.TGS_Gapcloser_Modified_V2_re.Gap_Closed/CuratedGapclosedASM.KOREF_S1.Trio_Specific.UpperCase.hap1.fasta",
        "hap2" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon/Final.TGS_Gapcloser_Modified_V2_re.Gap_Closed/CuratedGapclosedASM.KOREF_S1.Trio_Specific.UpperCase.hap2.fasta"
    }
    dict_ref_aft = {
        "hap1" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Reference_Fasta/Round1.SV_Polished_Fasta/PolishedASM.KOREF_S1.Trio_Specific.UpperCase.Polish_R1.SV.hap1.fasta",
        "hap2" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Reference_Fasta/Round1.SV_Polished_Fasta/PolishedASM.KOREF_S1.Trio_Specific.UpperCase.Polish_R1.SV.hap2.fasta"
    }

    path_output_format = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Reference_Fasta/Align_against_before/Align_Minimap2.{hap}.Ref_Bef.Qry_Aft.sorted.paf"

    for hap in ["hap1", "hap2"]:

        jobname = f"Align_BefAft_Polishing.Ref_Bef.Qry_Aft"
        
        run(
            dict_ref_aft[hap], 
            dict_ref_bef[hap], 
            path_output_format.format(hap = hap), 
            f"{dir_workspace}/Polishing/Minimap2_ASM_Align/tmp",
            f"{dir_workspace}/Polishing/Minimap2_ASM_Align/sh",
            f"{dir_workspace}/Polishing/Minimap2_ASM_Align/log",
            f"Minimap2.{jobname}.{curr_time}",
            48,
            [],
            [],
            **modified_option
        )