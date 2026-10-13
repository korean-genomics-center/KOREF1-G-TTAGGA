import os, time
import subprocess
from pathlib import Path

PATH_SCRIPT = os.path.join(str(Path(__file__).parents[0]), "check_alignment_of_matched_contig.py")


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

def run(path_paf, path_lowqv_region_info, path_output, anchor_min_idy, anchor_min_align, anchor_max_wronghang, anchor_min_align2hang, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):

    cmd = f"python {PATH_SCRIPT} --path_paf {path_paf} --path_lowqv_region_info {path_lowqv_region_info} --path_output {path_output} --anchor_min_idy {anchor_min_idy} --anchor_min_align {anchor_min_align} --anchor_max_wronghang {anchor_max_wronghang} --anchor_min_align2hang {anchor_min_align2hang} --n_jobs {n_jobs}"
    
    os.makedirs(os.path.dirname(path_output), exist_ok=True)
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
    
    curr_time = int(time.time())
    
    list_ref_assemblies = ['Hifiasm_Binned_HiFi', 'Hifiasm_Binned_ALT_HiFi', 'Hifiasm_Binned_HiFi_with0', 'Hifiasm_Binned_ALT_HiFi_with0', 'Flye_Binned_ONT', 'Flye_Binned_ONT.with0', 'Flye_Binned_HiFi_with0', 'Flye_Binned_HiFi']
    # list_ref_assemblies = ['Hifiasm_Binned_ALT_HiFi']

    dict_hap_to_parental = {
        1:'p',
        2:'m'
    }
    
    n_round = 5
    afterround = 2

    anchor_min_idy = 0.95
    anchor_min_align = 1000
    anchor_max_wronghang = 1000
    anchor_min_align2hang = 0.95
    
    for hap in [1, 2]:
        for name in list_ref_assemblies:
            path_paf = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Manual_Polishing/Round{n_round}_afterR{afterround}/Align.Binned_Contig/Align.REF_{name}.QRY_LowQVRegionHap{hap}_extra_range250.maphifiPw5s100cU101000000.filt_match1000_idy0.2.paf"
            path_lowqv_region_info = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Manual_Polishing/Round{n_round}_afterR{afterround}/Collect_LowQV_Region/LowQVRegion.KOREF_S1.Trio_Specific.UpperCase.Polish_R{n_round-1}.Man.hap{hap}.QV_per_window.k31.chunked.2000.qv_less60.extra_range250.region.tsv"

            path_output = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Manual_Polishing/Round{n_round}_afterR{afterround}/Anchored_Sequences/Filtered_Alignments/Hap{hap}/Filter_Anchor_Alignment.on.{name}.KOREF_S1.Trio_Specific.UpperCase.Polish_R{n_round-1}.Man.hap{hap}.extra_range250.idy95.a2h95.align1k.hang1k.tsv"

            assert os.path.exists(path_paf)
            assert os.path.exists(path_lowqv_region_info), path_lowqv_region_info
            run(
                path_paf, 
                path_lowqv_region_info, 
                path_output, 
                anchor_min_idy, 
                anchor_min_align, 
                anchor_max_wronghang, 
                anchor_min_align2hang, 
                f"{dir_workspace}/Polishing/Manual_Polishing/Filter_Anchor_Alignment/tmp", 
                f"{dir_workspace}/Polishing/Manual_Polishing/Filter_Anchor_Alignment/sh", 
                f"{dir_workspace}/Polishing/Manual_Polishing/Filter_Anchor_Alignment/log", 
                f"Filter_Anchor_Alignment.REF_{name}.QRY_KOREFPolishedASM.Round{n_round}_afterR{afterround}.Hap{hap}_Split_extra_range250.{curr_time}", 
                24, 
                [], 
                [], 
                **{}
            )
    
    