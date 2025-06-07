#%%
import json, os, sys, subprocess
from pathlib import Path

import pandas as pd
import numpy as np

sys.path.append(str(Path(__file__).parents[1]))
from fasta_parser import read_fasta_as_dict, save_dict_fasta_into_file

PATH_MINIMAP2 = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/minimap2/minimap2_2.29_r1283/minimap2-2.29_x64-linux/minimap2"
MINIMAP2_OPTION = "-x asm20 --secondary=no"

def run(path_ordered_json, path_fasta_input, dir_save, dir_workspace):
    dict_fasta_input = read_fasta_as_dict(path_fasta_input)
    dict_scaffold_ordered = read_json(path_ordered_json)

    for chrname in dict_scaffold_ordered.keys():
        dir_save_chrname = os.path.join(dir_save, chrname)
        os.makedirs(dir_save_chrname, exist_ok = True)

        list_scaffold_order = dict_scaffold_ordered[chrname]
        list_scaffold_order_nonna = list(filter(lambda val: val[0] != None, list_scaffold_order))

        for scfname, ori in list_scaffold_order_nonna:
            scf_sequence = dict_fasta_input[scfname].upper()
            if ori == '-':
                scf_sequence = get_reverse_complement_sequence(scf_sequence)
            save_dict_fasta_into_file({scfname:scf_sequence}, os.path.join(dir_save_chrname, f"{scfname}.fasta"))
        
        if len(list_scaffold_order_nonna) > 1:
            scfname_bef = list_scaffold_order_nonna[0][0]
            for scfname, _ in list_scaffold_order_nonna[1:]:
                align_option_suffix = "asm20_nosec"

                qryfasta =  os.path.join(dir_save_chrname, f"{scfname}.fasta")
                reffasta =  os.path.join(dir_save_chrname, f"{scfname_bef}.fasta")
                path_save = os.path.join(dir_save_chrname, "Align.Minimap2", f"Align.Ref.{scfname_bef}.QRY.{scfname}.{align_option_suffix}.paf")
                submit_minimap2_align(qryfasta, reffasta, path_save, f"{dir_workspace}/tmp", f"{dir_workspace}/sh", f"{dir_workspace}/log", f"Align.Ref.{scfname_bef}.QRY.{scfname}.{align_option_suffix}", 16)

def read_json(path_json):
    with open(path_json, 'r') as fr:
        dict_val = json.load(fr)
    return dict_val

def get_reverse_complement_sequence(sequence):
    dict_reverse = {
        'A' : 'T',
        'T' : 'A',
        'G' : 'C',
        'C' : 'G',
        'N' : 'N'
    }

    reversed_sequence = sequence[::-1]
    reverse_complement_sequence = list(map(dict_reverse.__getitem__, reversed_sequence))
    return ''.join(reverse_complement_sequence)

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

def submit_minimap2_align(qryfasta, reffasta, path_save, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):
    os.makedirs(os.path.dirname(path_save), exist_ok=True)

    cmd = f"{PATH_MINIMAP2} {MINIMAP2_OPTION} -o {path_save} -t {n_jobs}"
    
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

#%%
if __name__ == "__main__":
    hap = 2
    path_ordered_json = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling.RagTag.default_option.CHM13/CuratedDraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.non_t2t_scaffold.hap{hap}.RagTag.default_option.minimap2_latest.CHM13/CuratedDraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.hap{hap}.RagTag_with_Telomere_Assign.minimap2_latest.CHM13.pseudo_telomere_0.reoredered_scaffolds.json"
    path_fasta_input = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling/Adjusted/Telomere_Detection/Fasta_Split_by_Telomere/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg.no_telo.telo_single_masked.fa.k32.w100.z1000.ntLink.scaffolds.adjusted_telomere.no_telo.telo_single.fasta"
    dir_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling.RagTag.default_option.CHM13/CuratedDraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.non_t2t_scaffold.hap{hap}.RagTag.default_option.minimap2_latest.CHM13/Check_Overlap_of_Ordered_Scaffold"
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis/Minimap2_ASM_Align"

    run(path_ordered_json, path_fasta_input, dir_save, dir_workspace)
# %%
