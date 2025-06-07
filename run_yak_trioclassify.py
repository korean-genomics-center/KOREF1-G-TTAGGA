import pandas as pd
import os, time
import subprocess

PATH_YAK = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/yak/yak/yak"


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

def run(path_fastq, path_bin_pat, path_bin_mat, path_save_binned, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):
    cmd = f"""{PATH_YAK} triobin -t{n_jobs} {path_bin_pat} {path_bin_mat} {path_fastq} > {path_save_binned}"""
    
    os.makedirs(os.path.dirname(path_save_binned), exist_ok=True)
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

    path_fastq = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling/Adjusted/Parental_Split/Check_classify/h1tg000019l_tail_74975.fasta"
    paternal_count = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/PATERNAL/Yak_Binning/PATERNAL_Trimmed_F15_T5.noNbase.merged.yak"
    maternal_count = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/MATERNAL/Yak_Binning/MATERNAL_Trimmed_F15_T5.noNbase.merged.yak"
    path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling/Adjusted/Parental_Split/Check_classify/h1tg000019l_tail_74975.trio_classify.txt"
    
    curr_time = int(time.time())
    run(
        path_fastq,
        paternal_count,
        maternal_count,
        path_save,
        f"{dir_workspace}/Classify_Scaffold_by_Trio/tmp",
        f"{dir_workspace}/Classify_Scaffold_by_Trio/sh",
        f"{dir_workspace}/Classify_Scaffold_by_Trio/log",
        f"Classify_Scaffold_by_Trio.KOREF.KOREF_S1.ul_L100kb_QV15.Trio_Specific.h1tg000019l_tail_74975.{curr_time}",
        16,
        [],
        []
    ) 


if __name__ == "__main2__":
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis"
    
    
    curr_time = int(time.time())
    
    for hap in [1, 2]:
        path_fastq = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling/Adjusted/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg.no_telo.telo_single_masked.fa.k32.w100.z1000.ntLink.scaffolds.adjusted_telomere.fa"
        paternal_count = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/PATERNAL/Yak_Binning/PATERNAL_Trimmed_F15_T5.noNbase.merged.yak"
        maternal_count = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/MATERNAL/Yak_Binning/MATERNAL_Trimmed_F15_T5.noNbase.merged.yak"
        
        path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling/Adjusted/Parental_Split/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg.no_telo.telo_single_masked.fa.k32.w100.z1000.ntLink.scaffolds.adjusted_telomere.split_info.parental_specific_kmer.txt"
        
        curr_time = int(time.time())
        run(
            path_fastq,
            paternal_count,
            maternal_count,
            path_save,
            f"{dir_workspace}/Classify_Scaffold_by_Trio/tmp",
            f"{dir_workspace}/Classify_Scaffold_by_Trio/sh",
            f"{dir_workspace}/Classify_Scaffold_by_Trio/log",
            f"Classify_Scaffold_by_Trio.KOREF.KOREF_S1.ul_L100kb_QV15.Trio_Specific.hap{hap}.ntLink.adjust_telomere.Parental_Specifric.{curr_time}",
            16,
            [],
            []
        )

if __name__ == "__main2__":
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis"
    
    
    curr_time = int(time.time())
    
    for hap in [1, 2]:
        path_fastq = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling/Adjusted/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.asm.dip.hap{hap}.p_ctg.no_telo.telo_single_masked.fa.k32.w100.z1000.ntLink.scaffolds.adjusted_telomere.fa"
        paternal_count = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/PATERNAL/Yak_Binning/PATERNAL_Trimmed_F15_T5.noNbase.merged.yak"
        maternal_count = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/WGS/Trio/MATERNAL/Yak_Binning/MATERNAL_Trimmed_F15_T5.noNbase.merged.yak"
        
        path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling/Adjusted/Parental_Split/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.asm.dip.hap{hap}.p_ctg.no_telo.telo_single_masked.fa.k32.w100.z1000.ntLink.scaffolds.adjusted_telomere.split_info.parental_specific_kmer.txt"
        
        curr_time = int(time.time())
        run(
            path_fastq,
            paternal_count,
            maternal_count,
            path_save,
            f"{dir_workspace}/Classify_Scaffold_by_Trio/tmp",
            f"{dir_workspace}/Classify_Scaffold_by_Trio/sh",
            f"{dir_workspace}/Classify_Scaffold_by_Trio/log",
            f"Classify_Scaffold_by_Trio.KOREF.KOREF_S1.ul_L100kb_QV15.hap{hap}.ntLink.adjust_telomere.Parental_Specifric.{curr_time}",
            16,
            [],
            []
        )