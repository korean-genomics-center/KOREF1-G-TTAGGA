#%%
import subprocess, os, time
from pathlib import Path

PATH_SCRIPT = os.path.join(os.path.dirname(__file__), "get_qv_only.sh")
DIR_MERQURY = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/merqury"
DIR_MERYL = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/meryl-1.4.1/bin"
DIR_SAMTOOLS = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/Samtools/samtools-1.17/samtools-1.17/bin"
DIR_BEDTOOLS = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/Bedtools"

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

def make_symlink(path_src, path_dst):
    if os.path.exists(path_dst):
        # return
        os.remove(path_dst)
    os.symlink(path_src, path_dst)

def run(path_asm, path_child_meryl, path_output_prefix, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):
    path_dir_save = os.path.dirname(path_output_prefix)
    os.makedirs(path_dir_save, exist_ok=True)
    output_prefix = os.path.basename(path_output_prefix)

    path_asm_link = os.path.join(path_dir_save, os.path.basename(path_asm))
    make_symlink(path_asm, path_asm_link)

    cmd = f"""cd {path_dir_save};{PATH_SCRIPT} {path_child_meryl} {path_asm_link} {output_prefix}"""
    
    merqury_environ = os.environ.copy()
    merqury_environ["PATH"] = f"{DIR_MERYL}:{DIR_SAMTOOLS}:{DIR_BEDTOOLS}:{merqury_environ['PATH']}"
    merqury_environ["MERQURY"] = DIR_MERQURY
    merqury_environ["OMP_NUM_THREADS"] = f"{n_jobs}"
    
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
        ),
        env = merqury_environ
    )

if __name__ == "__main2__":
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis"
    
    dict_path_assembly_format = {
        # Hifiasm Binned (HiFi)
        "Hifiasm_Binned_HiFi":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20/Hap{hap}/KOREF_S1.Binned_Assembly.hifi_L1000_QV20.triobin_{triobin}.haploid.l0.asm.bp.p_ctg.fasta",

        # Hifiasm Binned (Alt) (HiFi)
        "Hifiasm_Binned_ALT_HiFi":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20/Hap{hap}/KOREF_S1.Binned_Assembly.hifi_L1000_QV20.triobin_{triobin}.haploid.l0.asm.bp.a_ctg.fasta",

        # Hifiasm Binned with 0 (HiFi)
        "Hifiasm_Binned_with0_HiFi":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20.with_0/Hap{hap}/KOREF_S1.Binned_Assembly.hifi_L1000_QV20.triobin_{triobin}_0.haploid.l0.asm.bp.p_ctg.fasta",

        # Hifiasm Binned with 0 (Alt) (HiFi)
        "Hifiasm_Binned_with0_ALT_HiFi":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20.with_0/Hap{hap}/KOREF_S1.Binned_Assembly.hifi_L1000_QV20.triobin_{triobin}_0.haploid.l0.asm.bp.a_ctg.fasta",

        # Flye Binned (ONT)
        "Flye_Binned_ONT":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/ONT_UL.Filter_L50kb_QV20_Q30Perc80/Hap{hap}/assembly.fasta",
        
        # # Flye Binned with 0 (ONT)
        # "Flye_Binned_ONT.with0":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/ONT_UL.Filter_L50kb_QV20_Q30Perc80.with_0/Hap{hap}/assembly.fasta",

        # Flye Binned (HiFi)
        "Flye_Binned_HiFi":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20/Hap{hap}/assembly.fasta",

        # Flye Binned with 0 (HiFi)
        "Flye_Binned_HiFi.with0":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20.with_0/Hap{hap}/assembly.fasta",

    #     # "Unplaced_Scaffold" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/Chromosome_level_scaffolds.Assign_by_CHM13/CuratedDraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.hap{hap}.RagTag_with_Telomere_Assign.minimap2_latest.CHM13.unplaced.fasta"

    }

    dict_path_assembly_format = {
        "Polished_ASM":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Reference_Fasta/Round8.SV_SNV_Polished_Fasta/PolishedASM.KOREF_S1.Trio_Specific.UpperCase.Polish_R8.SV_SNV.hap{hap}.fasta"
    }

    k = 21
    
    path_child_meryl = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/PacBio/Meryl_Binning/PacBio_HiFi_ccs_converted.filt.Filter_L1000_QV20.k{k}.meryl"
    
    dict_hap_to_parental = {
        1:'p',
        2:'m'
    }

    curr_time = int(time.time())

    for hap in [1, 2]:
        for asm_name, path_asm_format in dict_path_assembly_format.items():
            path_asm = path_asm_format.format(hap = hap, triobin = dict_hap_to_parental[hap])
            # path_asm = path_asm_format

            assert os.path.exists(path_asm), path_asm

            dir_asm = os.path.dirname(path_asm)
            filename_asm = os.path.basename(path_asm)

            outputname_prefix = filename_asm.replace(".fasta", f".k{k}")

            path_output_prefix = os.path.join(dir_asm, "Merqury_QV", outputname_prefix)

            run(
                path_asm,
                path_child_meryl,
                path_output_prefix,
                f"{dir_workspace}/Merqury_Evaluation/tmp",
                f"{dir_workspace}/Merqury_Evaluation/sh",
                f"{dir_workspace}/Merqury_Evaluation/log",
                f"Evaluate_QV_only.KOREF_S1.{asm_name}.Hap{hap}.k{k}.{curr_time}",
                # f"Evaluate_QV_only.KOREF_S1.{asm_name}.k21.{curr_time}",
                30,
                [],
                []
            )
#%%
if __name__ == "__main__":
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis"
    
    k = 21
    
    curr_time = int(time.time())
    
    path_child_meryl = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/PacBio/Meryl_Binning/PacBio_HiFi_ccs_converted.filt.Filter_L1000_QV20.k{k}.meryl"
    
    path_asm = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Final_Assembly/FinalASM.KOREF_S1.Trio_Specific.Combined.Reordered.fasta"

    path_output_prefix = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Merqury_Evaluation/Merqury_QV/FinalASM.KOREF_S1.Trio_Specific.Combined.Reordered.k{k}"

    run(
        path_asm,
        path_child_meryl,
        path_output_prefix,
        f"{dir_workspace}/Merqury_Evaluation/tmp",
        f"{dir_workspace}/Merqury_Evaluation/sh",
        f"{dir_workspace}/Merqury_Evaluation/log",
        f"Evaluate_QV_only.KOREF_S1.Combined.k{k}.{curr_time}",
        # f"Evaluate_QV_only.KOREF_S1.{asm_name}.k21.{curr_time}",
        30,
        [],
        []
    )