import pandas as pd
import os, time
import subprocess

PATH_GFASTATS = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/gfastats"


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

def run_gfastats(path_assembly, path_save, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):
    cmd = f"""{PATH_GFASTATS} {path_assembly} -t -j {n_jobs} > {path_save}"""
    
    for param_name, param_value in additional_params.items():
        cmd += f" --{param_name} {param_value}"
    
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



if __name__ == "__main2__":
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis/Stat_Assembly"
    
    
    dict_hap_to_parental = {
        1:'p',
        2:'m'
    }
    
    curr_time = int(time.time())
    for hap in [1, 2]:
    
        path_assembly = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20.with_0/Hap{hap}/KOREF_S1.Binned_Assembly.hifi_L1000_QV20.triobin_{dict_hap_to_parental[hap]}_0.haploid.l0.asm.bp.p_ctg.fasta"
        path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20.with_0/Hap{hap}/KOREF_S1.Binned_Assembly.hifi_L1000_QV20.triobin_{dict_hap_to_parental[hap]}_0.haploid.l0.asm.bp.p_ctg.fasta.stats"
        
        curr_time = int(time.time())
        run_gfastats(
            path_assembly,
            path_save,
            f"{dir_workspace}/Stat_Assembly/tmp",
            f"{dir_workspace}/Stat_Assembly/sh",
            f"{dir_workspace}/Stat_Assembly/log",
            f"Stat_Assembly.Hifiasm.KOREF_S1.HiFi.Binned_Assembly_with0.hap{hap}.{curr_time}",
            16,
            [],
            []
        )


if __name__ == "__main__":
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis/Stat_Assembly"
    
    dict_path_assembly_format = {
        # Hifiasm Binned (HiFi)
        "Hifiasm_Binned_HiFi":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20/Hap{hap}/KOREF_S1.Binned_Assembly.hifi_L1000_QV20.triobin_{triobin}.haploid.l0.asm.bp.p_ctg.qv50.hamming0.05.fasta",

        # Hifiasm Binned (Alt) (HiFi)
        "Hifiasm_Binned_ALT_HiFi":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20/Hap{hap}/KOREF_S1.Binned_Assembly.hifi_L1000_QV20.triobin_{triobin}.haploid.l0.asm.bp.a_ctg.qv50.hamming0.05.fasta",

        # Hifiasm Binned with 0 (HiFi)
        "Hifiasm_Binned_with0_HiFi":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20.with_0/Hap{hap}/KOREF_S1.Binned_Assembly.hifi_L1000_QV20.triobin_{triobin}_0.haploid.l0.asm.bp.p_ctg.qv50.hamming0.05.fasta",

        # Hifiasm Binned with 0 (Alt) (HiFi)
        "Hifiasm_Binned_with0_ALT_HiFi":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20.with_0/Hap{hap}/KOREF_S1.Binned_Assembly.hifi_L1000_QV20.triobin_{triobin}_0.haploid.l0.asm.bp.a_ctg.qv50.hamming0.05.fasta",

        # # Flye Binned (ONT)
        # "Flye_Binned_ONT":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/ONT_UL.Filter_L50kb_QV20_Q30Perc80/Hap{hap}/assembly.fasta",

        # Flye Binned (HiFi)
        "Flye_Binned_HiFi_with0":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20/Hap{hap}/assembly.qv50.hamming0.05.fasta",

        # Flye Binned with 0 (HiFi)
        "Flye_Binned_HiFi":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20.with_0/Hap{hap}/assembly.qv50.hamming0.05.fasta",

        # "Unplaced_Scaffold" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/Chromosome_level_scaffolds.Assign_by_CHM13/CuratedDraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.hap{hap}.RagTag_with_Telomere_Assign.minimap2_latest.CHM13.unplaced.fasta"

    }
    
    dict_hap_to_parental = {
        1:'p',
        2:'m'
    }
    
    curr_time = int(time.time())
    for hap in [1, 2]:
        for asm_name, path_asm_format in dict_path_assembly_format.items():
            path_asm = path_asm_format.format(hap = hap, triobin = dict_hap_to_parental[hap])
    
            path_save = path_asm.replace(".fasta", ".fasta.stats")
            assert os.path.exists(path_asm), path_asm
            curr_time = int(time.time())
            run_gfastats(
                path_asm,
                path_save,
                f"{dir_workspace}/Stat_Assembly/tmp",
                f"{dir_workspace}/Stat_Assembly/sh",
                f"{dir_workspace}/Stat_Assembly/log",
                f"Stat_Assembly.{asm_name}.hap{hap}.qv50hamming0.05.{curr_time}",
                16,
                [],
                []
            )