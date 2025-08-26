import os, time
import subprocess

PATH_MINIMAP2 = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/minimap2/minimap2_2.29_r1283/minimap2-2.29_x64-linux/minimap2"


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

def run(qryfasta, reffasta, path_save, min_match_length, minimap2_option, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):
    os.makedirs(os.path.dirname(path_save), exist_ok=True)

    cmd = f"{PATH_MINIMAP2} {minimap2_option} -t {n_jobs}"
    
    for param_name, param_value in additional_params.items():
        if len(param_name) == 1:
            cmd += f" -{param_name} {param_value}"
        else:
            cmd += f" --{param_name} {param_value}"
    cmd += f" {reffasta} {qryfasta} | awk '{{if($10>{min_match_length}) print $0}}' - > {path_save}"
    
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
    
    dict_binnedasm_format = {
        # Hifiasm Binned (HiFi)
        "Hifiasm_Binned_HiFi":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20/Hap{hap}/KOREF_S1.Binned_Assembly.hifi_L1000_QV20.triobin_{triobin}.haploid.l0.asm.bp.p_ctg.fasta",

        #Hifiasm Binned (Alt) (HiFi)
        "Hifiasm_Binned_ALT_HiFi":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20/Hap{hap}/KOREF_S1.Binned_Assembly.hifi_L1000_QV20.triobin_{triobin}.haploid.l0.asm.bp.a_ctg.fasta",

        # Hifiasm Binned with0 (HiFi)
        "Hifiasm_Binned_HiFi_with0":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20.with_0/Hap{hap}/KOREF_S1.Binned_Assembly.hifi_L1000_QV20.triobin_{triobin}_0.haploid.l0.asm.bp.p_ctg.fasta",
        
        #Hifiasm Binned with0 (Alt) (HiFi)
        "Hifiasm_Binned_ALT_HiFi_with0":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20.with_0/Hap{hap}/KOREF_S1.Binned_Assembly.hifi_L1000_QV20.triobin_{triobin}_0.haploid.l0.asm.bp.a_ctg.fasta",

        # Flye Binned with0 (HiFi)
        "Flye_Binned_HiFi_with0":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20.with_0/Hap{hap}/assembly.fasta",

        # Flye Binned (HiFi)
        "Flye_Binned_HiFi":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20/Hap{hap}/assembly.fasta",

        # Unplaced Scaffold
        "Unplaced_Scaffold" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/Chromosome_level_scaffolds.Assign_by_CHM13/CuratedDraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.hap{hap}.RagTag_with_Telomere_Assign.minimap2_latest.CHM13.unplaced.fasta"
    }

    dict_hap_to_parental = {
        1:'p',
        2:'m'
    }

    min_match_length = 100
    minimap2_option = "-x asm5 -P -c -K 5G"
    
    for hap in [1, 2]:
        triobin = dict_hap_to_parental[hap]
        
        list_binnedasm = list(dict_binnedasm_format.keys())
        for name1 in list_binnedasm:
            for name2 in list_binnedasm:
            
                path_fasta1 = dict_binnedasm_format[name1].replace("{hap}", str(hap)).replace("{triobin}", triobin)
                path_fasta2 = dict_binnedasm_format[name2].replace("{hap}", str(hap)).replace("{triobin}", triobin)
                
                path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon_wo_RawONT/Contig_Alignment/Hap{hap}/Align.Pairwise/Align.Hap{hap}.REF.{name1}.QRY.{name2}.asm5Pc_triobin_{triobin}.filt_match{min_match_length}.paf"

                assert os.path.exists(path_fasta1), path_fasta1
                assert os.path.exists(path_fasta2), path_fasta2
                run(
                    path_fasta2, 
                    path_fasta1, 
                    path_save, 
                    min_match_length,
                    minimap2_option,
                    f"{dir_workspace}/Minimap2_ASM_Align/tmp",
                    f"{dir_workspace}/Minimap2_ASM_Align/sh",
                    f"{dir_workspace}/Minimap2_ASM_Align/log",
                    f"Minimap2.Hap{hap}.REF_{name1}_{triobin}.QRY_{name2}_{triobin}.asm5Pc.{curr_time}",
                    16,
                    [],
                    [],
                    **{}
                )

    