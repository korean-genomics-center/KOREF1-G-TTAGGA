import pandas as pd
import os, time
import subprocess

PATH_TGSSEQSPLIT = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/TGSGapcloser/TGSGapcloser_1.2.1/TGS-GapCloser/tgsgapcloserbin/tgsseqsplit"


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

def run(path_scaffolds, path_output_prefix, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list()):    

    cmd = f"""{PATH_TGSSEQSPLIT} --input_scaff {path_scaffolds} --prefix {path_output_prefix}"""
        
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
    
    for hap in [1, 2]:
        for telorepeat in [5]:
            path_scaffolds = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/Round4.TGS_Gapcloser_Modified.Unplaced_Scf_asm5P_idy8_match1000_hard.HiFi_L1000_QV20_parental_maphifiN100P_hard.Hifiasm_BinnedASM_asm5P_hard.Flye_ONT_UL_BinnedASM_asm5P_hard/Hap{hap}/n_pseudo_telorepeat{telorepeat}/CuratedDraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.hap{hap}.Assigned_with_CHM13.n_pseudo_telorepeat{telorepeat}.GapFilling.Unplaced_Scf.HiFi_L1000_QV20_parental_{dict_hap_to_parental[hap]}.Hifiasm_BinnedASM.Flye_ONT_UL_BinnedASM.scaff_seqs"

            path_output_prefix = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/Round5.TGS_Gapcloser_Modified.Unplaced_Scf_asm5P_idy8_match1000_hard.HiFi_L1000_QV20_parental_maphifiN100P_hard.Hifiasm_BinnedASM_asm5P_hard.Flye_ONT_UL_BinnedASM_asm5P_hard.Manual_Scaffolding/Prealignment_for_Candidate_Filtering/Hap{hap}/CuratedDraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.hap{hap}.Assigned_with_CHM13.n_pseudo_telorepeat{telorepeat}.GapFilling.Unplaced_Scf.HiFi_L1000_QV20_parental_{dict_hap_to_parental[hap]}.Hifiasm_BinnedASM.Flye_ONT_UL_BinnedASM"

            curr_time = int(time.time())
            run(path_scaffolds, 
                path_output_prefix, 
                f"{dir_workspace}/Split_Scaffold/tmp", 
                f"{dir_workspace}/Split_Scaffold/sh", 
                f"{dir_workspace}/Split_Scaffold/log", 
                f"Split_Scaffold.KOREF_S1.ul_L100kb_QV15.Trio_Specific.hap{hap}.n_pseudo_telorepeat{telorepeat}.R5.Manual_Scaffolding.asm5.{curr_time}", 
                1
            )
#%%