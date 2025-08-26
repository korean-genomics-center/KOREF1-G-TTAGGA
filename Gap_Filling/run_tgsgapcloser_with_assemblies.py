import pandas as pd
import os, time
import subprocess

PATH_TGSGAPCLOSER = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/TGSGapcloser/TGSGapcloser_1.2.1_Modifiy/TGS-GapCloser/tgsgapcloser"


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

def run(path_scaffolds, path_reads, path_output_prefix, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list()):    

    cmd = f"""cd {os.path.dirname(path_output_prefix)};{PATH_TGSGAPCLOSER} --scaff {path_scaffolds} --reads {path_reads} --output {path_output_prefix} --ne --minmap_arg ' -x asm5 -P -c -K 4G' --tgstype pb --min_idy 0.9 --min_match 500 --hard --min_idy_hard 0.9 --min_align_hard 500 --max_hang_hard 2000  --min_alignratio_hard 0.9 --thread {n_jobs}"""
        
    os.makedirs(os.path.dirname(path_output_prefix), exist_ok=True)
    os.makedirs(dir_tmp, exist_ok=True)
    os.makedirs(dir_shell, exist_ok=True)
    os.makedirs(dir_log, exist_ok=True)
    
    # my_env = os.environ.copy()
    # my_env["PATH"] = f"{PATH_MINIMAP2}:{my_env['PATH']}"
    
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
        env = None
    )



if __name__ == "__main__":
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis"

    dict_hap_to_parental = {
        1:'p',
        2:'m'
    }

    dict_hap_to_local = {
        1:"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/Local_Assembly/HiFi/Hap1/chr13_1_5000000.Mapped_KOREF_DraftASM_hap2.with_0.l0/KOREF_S1.Local_Assembly.chr13_1_5000000.hifi_L1000_QV20.triobin_p_0.haploid.l0.asm.bp.p_ctg.fasta",
        2:"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/Local_Assembly/HiFi/Hap2/chr15_94775221.Mapped_KOREF_DraftASM_hap1.with_0.l0/KOREF_S1.Local_Assembly.chr15_94775221.hifi_L1000_QV20.triobin_m_0.haploid.l0.asm.bp.p_ctg.fasta"
    }
    
    for hap in [1]:
        path_scaffolds = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon_wo_RawONT/Round7.Flye_Hifi_Bin0/Hap{hap}/CuratedDraftASM.KOREF_S1.Trio_Specific.hap{hap}.Round7.GapFilling.Manual_Scaff.LocalASM.Unplaced_Scf.Hifiasm_Hifi_Bin.Hifiasm_Hifi_Bin0.Flye_Hifi_Bin.Flye_Hifi_Bin0.scaff_seqs"

        # Local assembly
        # path_reads = dict_hap_to_local[hap]

        # Hifiasm Binned (HiFi)
        # path_reads = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20/Hap{hap}/KOREF_S1.Binned_Assembly.hifi_L1000_QV20.triobin_{dict_hap_to_parental[hap]}.haploid.l0.asm.bp.p_ctg.qv50.hamming0.05.fasta"

        # Hifiasm Binned with 0 (HiFi)
        # path_reads = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20.with_0/Hap{hap}/KOREF_S1.Binned_Assembly.hifi_L1000_QV20.triobin_{dict_hap_to_parental[hap]}_0.haploid.l0.asm.bp.p_ctg.qv50.hamming0.05.fasta"

        # Flye Binned (HiFi)
        # path_reads = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20/Hap{hap}/assembly.qv50.hamming0.05.fasta"

        # Flye Binned with 0 (HiFi)
        # path_reads = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20.with_0/Hap{hap}/assembly.qv50.hamming0.05.fasta"


        # Unplaced Scaffolds
        # path_reads = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/Chromosome_level_scaffolds.Assign_by_CHM13/CuratedDraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.hap{hap}.RagTag_with_Telomere_Assign.minimap2_latest.CHM13.unplaced.fasta"
        

        path_output_prefix = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon_wo_RawONT/Round8.Hifiasm_Hifi_Bin0_a/Hap{hap}/CuratedDraftASM.KOREF_S1.Trio_Specific.hap{hap}.Round8.GapFilling.Manual_Scaff.LocalASM.Unplaced_Scf.Hifiasm_Hifi_Bin.Hifiasm_Hifi_Bin0.Flye_Hifi_Bin.Flye_Hifi_Bin0.Hifiasm_Hifi_Bin0_a"

        assert os.path.exists(path_scaffolds), path_scaffolds
        assert os.path.exists(path_reads), path_reads
        
        curr_time = int(time.time())
        run(path_scaffolds,                                 
            path_reads, 
            path_output_prefix, 
            f"{dir_workspace}/GapFilling/tmp", 
            f"{dir_workspace}/GapFilling/sh", 
            f"{dir_workspace}/GapFilling/log", 
            f"Gapfilling_Mod_V2_re.KOREF_S1.ul_L100kb_QV15.Trio_Specific.hap{hap}.R8.Hifiasm_Hifi_Bin0_a.{curr_time}", 
            96
        )
#%%