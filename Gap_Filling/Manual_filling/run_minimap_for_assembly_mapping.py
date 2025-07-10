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

def run(qryfasta, reffasta, path_save, min_match, min_idy, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):
    os.makedirs(os.path.dirname(path_save), exist_ok=True)

    cmd = f"{PATH_MINIMAP2} -x asm20 -k15 -w5 -e0 -U1,1000000 -P --no-long-join -r500,20k --q-occ-frac 0 -f 0.0001 -K 150M -t {n_jobs} "
    # cmd = f"{PATH_MINIMAP2} -x asm5 -P -c -K 3G -o {path_save} -t {n_jobs}  -s100 -z800,100 "
    
    for param_name, param_value in additional_params.items():
        if len(param_name) == 1:
            cmd += f" -{param_name} {param_value}"
        else:
            cmd += f" --{param_name} {param_value}"
    cmd += f" {reffasta} {qryfasta}"
    cmd += f"| awk '{{if(($10>{min_match}) && (($10/$11)>{min_idy})) print $0}}' - > {path_save}"
    
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
            int(n_jobs*1.5),
            hostname=hostname,
            hold_jids=jobs_wait
        )
    )

if __name__ == "__main__":
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis"
    
    curr_time = int(time.time())
    
    dict_binnedasm_format = {
        "Hifiasm_Binned_pri":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.Binned_Input/ONT_UL.Filter_L50kb_QV20_Q30Perc80.HiFi.Filter_L1000_QV20/Hap{hap}/KOREF_S1.Binned_Assembly.ul_L50kb_QV20_Q30Perc80.triobin_{triobin}.hifi_L1000_QV20.triobin_{triobin}.haploid.l0.asm.bp.p_ctg.fasta",
        "Hifiasm_Binned_alt":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.Binned_Input/ONT_UL.Filter_L50kb_QV20_Q30Perc80.HiFi.Filter_L1000_QV20/Hap{hap}/KOREF_S1.Binned_Assembly.ul_L50kb_QV20_Q30Perc80.triobin_{triobin}.hifi_L1000_QV20.triobin_{triobin}.haploid.l0.asm.bp.a_ctg.fasta",
        "Flye_HiFi_Binned":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20/Hap{hap}/assembly.fasta",
        "Flye_ONT_Binned" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/ONT_UL.Filter_L50kb_QV20_Q30Perc80/Hap{hap}/assembly.fasta",
        "Unplaced_Scaffold" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/Chromosome_level_scaffolds.Assign_by_CHM13/CuratedDraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.hap{hap}.RagTag_with_Telomere_Assign.minimap2_latest.CHM13.unplaced.fasta",
        "Flye_HiFi_Binned_with0" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20.with_0/Hap{hap}/assembly.fasta",
        "Flye_ONT_Binned_with0" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/ONT_UL.Filter_L50kb_QV20_Q30Perc80.with_0/Hap{hap}/assembly.fasta",
        "Flye_HiFi_Binned_with0a" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20.with_0a/Hap{hap}/assembly.fasta",
        "Flye_ONT_Binned_with0a" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/ONT_UL.Filter_L50kb_QV20_Q30Perc80.with_0a/Hap{hap}/assembly.fasta",
        "Flye_ONT_Binned_short" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/ONT_UL.Filter_L10kb_QV15/Hap{hap}/assembly.fasta"
    }

    dict_hap_to_parental = {
        1:'p',
        2:'m'
    }
    
    min_match = 100
    min_idy = 0.2
    
    for hap in [1, 2]:
        for name, path_format in dict_binnedasm_format.items():
            if name != "Flye_HiFi_Binned":
                continue
            triobin = dict_hap_to_parental[hap]
            
            path_ref = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon/Round13.TGS_Gapcloser_Modified_V2_re.Unpl_Scf.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.Manual_Scaff.HiFi_trio.Hifiasm_Loc.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.ONT_trio.ONT_Amb.Manual_Scaff/Prealignment_for_Candidate_Filtering/Hap{hap}/CuratedDraftASM.KOREF_S1.Trio_Specific.hap{hap}.GapFilling.Unplaced_Scf.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.Manual_Scaff.HiFi_trio.Hifiasm_Loc.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.ONT_trio.ONT_Amb.Manual_Scaff.contig"

            path_fasta = path_format.replace("{hap}", str(hap)).replace("{triobin}", triobin)
            
            path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon/Round13.TGS_Gapcloser_Modified_V2_re.Unpl_Scf.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.Manual_Scaff.HiFi_trio.Hifiasm_Loc.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.ONT_trio.ONT_Amb.Manual_Scaff/Prealignment_for_Candidate_Filtering/Hap{hap}/Align.Option_Adjust/Align.REF_KOREFDraftASM_Hap{hap}_Split.QRY_{name}.asm20Pk15w5e0U1r50020knolongjoinqoccfrac0f0001.filt_match{min_match}_idy{min_idy}.paf"

            assert os.path.exists(path_fasta)
            assert os.path.exists(path_ref)
            run(
                path_fasta, 
                path_ref, 
                path_save, 
                min_match,
                min_idy,
                f"{dir_workspace}/Minimap2_ASM_Align/tmp",
                f"{dir_workspace}/Minimap2_ASM_Align/sh",
                f"{dir_workspace}/Minimap2_ASM_Align/log",
                f"Minimap2.REF_KOREFDraftASM_Hap{hap}_Split.QRY_{name}.asm20Pk15w5e0U1r50020knolongjoinqoccfrac0f0001.{curr_time}",
                24,
                [],
                [],
                **{}
            )
    
    