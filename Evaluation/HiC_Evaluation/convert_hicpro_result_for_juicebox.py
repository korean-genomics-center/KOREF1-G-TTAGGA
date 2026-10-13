#%%
import subprocess, os, time

PATH_HICPRO2JUICEBOX = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/HiC-Pro/HiC-Pro_v3.1.0/HiC-Pro/bin/utils/hicpro2juicebox.sh"
PATH_JUICERTOOLS = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/Juicer/Juicer_v1.6_master/juicer_1.6_CPU/juicer_tools_1.22.01.jar"


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

def write_shell_script_and_run(path_shell, shell_context, env = None):
    with open(path_shell, 'w') as fw:
        fw.write(shell_context)
    
    subprocess.run(f"qsub {path_shell}", shell = True, env = env)

def run(path_validpairs, path_genomesize, path_fragment_bed, dir_output, dir_save_tmp, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):
    
    cmd = f"""{PATH_HICPRO2JUICEBOX} -i {path_validpairs} -g {path_genomesize} -j {PATH_JUICERTOOLS} -r {path_fragment_bed} -o {dir_output} -t {dir_save_tmp} -e {n_jobs}"""
    
    os.makedirs(dir_output, exist_ok=True)
    os.makedirs(dir_save_tmp, exist_ok=True)
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
    
    
    for hap in [1, 2]:
        path_validpairs = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/HiC_Evaluation/HiC_Alignment/NovaSeq_20190425_DpnII-HinfI.MapQ0/Hap{hap}/results/hic_results/data/sample/KOREF-S1_PBMC_Hi-C_NovaSeq_20190425_DpnII-HinfI_L01.dedup.trim_F15_T0_FinalASM.KOREF_S1.Trio_Specific.hap{hap}.Reordered.bwt2pairs.validPairs"
        path_genomesize = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/HiC_Evaluation/HiC_Alignment/NovaSeq_20190425_DpnII-HinfI.MapQ0/Hap{hap}/genome/FinalASM.KOREF_S1.Trio_Specific.hap{hap}.Reordered.fasta.genome"
        path_fragment_bed = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/HiC_Evaluation/HiC_Alignment/NovaSeq_20190425_DpnII-HinfI.MapQ0/Hap{hap}/genome/FinalASM.KOREF_S1.Trio_Specific.hap{hap}.Reordered.fasta.GATC_GANTC.bed"

        dir_save_tmp = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/HiC_Evaluation/HiC_Alignment/NovaSeq_20190425_DpnII-HinfI.MapQ0/Hap{hap}/HiCPro2Juicer/tmp"
        dir_output = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/HiC_Evaluation/HiC_Alignment/NovaSeq_20190425_DpnII-HinfI.MapQ0/Hap{hap}/HiCPro2Juicer/"
        
        curr_time = int(time.time())
        run(
            path_validpairs, path_genomesize, path_fragment_bed, dir_output, dir_save_tmp,
            f"{dir_workspace}/HiC_Processing/tmp",
            f"{dir_workspace}/HiC_Processing/sh",
            f"{dir_workspace}/HiC_Processing/log",
            f"HiCPro2Juicebox.KOREF_S1.FinalASM.Hap{hap}.{curr_time}",
            n_jobs = 64,
            hostname = [],
            jobs_wait = []
        )
#%%