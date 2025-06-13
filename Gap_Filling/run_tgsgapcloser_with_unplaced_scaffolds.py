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

    cmd = f"""cd {os.path.dirname(path_output_prefix)};{PATH_TGSGAPCLOSER} --scaff {path_scaffolds} --reads {path_reads} --output {path_output_prefix} --ne --minmap_arg ' -x asm5 -P' --tgstype pb --min_idy 0.8 --min_match 1000 --hard --min_idy_hard 0.8 --min_match_hard 2000 --max_wronghang_hard 2000 --thread {n_jobs}"""
        
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

    for hap in [1, 2]:
        for telorepeat in [5]:
            path_scaffolds = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/Chromosome_level_scaffolds.Assign_by_CHM13/CuratedDraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.with_teloboth.hap{hap}.RagTag_with_Telomere_Assign.minimap2_latest.CHM13.n_pseudo_telorepeat{telorepeat}.fasta"
            path_reads = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/Chromosome_level_scaffolds.Assign_by_CHM13/CuratedDraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.hap{hap}.RagTag_with_Telomere_Assign.minimap2_latest.CHM13.unplaced.fasta"

            path_output_prefix = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/Round1.TGS_Gapcloser_Modified.Unplaced_Scf_asm5P_idy8_match1000_hard/Hap{hap}/n_pseudo_telorepeat{telorepeat}/CuratedDraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.hap{hap}.Assigned_with_CHM13.n_pseudo_telorepeat{telorepeat}.GapFilling.Unplaced_Scf"

            curr_time = int(time.time())
            run(path_scaffolds, 
                path_reads, 
                path_output_prefix, 
                f"{dir_workspace}/GapFilling/tmp", 
                f"{dir_workspace}/GapFilling/sh", 
                f"{dir_workspace}/GapFilling/log", 
                f"Gapfilling_Mod.KOREF_S1.ul_L100kb_QV15.Trio_Specific.hap{hap}.n_pseudo_telorepeat{telorepeat}.R1.Unplaced_Scf.asm5P_hard.{curr_time}", 
                16
            )
#%%