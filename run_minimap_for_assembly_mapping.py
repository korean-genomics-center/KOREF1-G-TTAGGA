import os, time
import subprocess

PATH_MINIMAP2 = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/minimap2/minimap2_2.29_r1283/minimap2-2.29_x64-linux/minimap2"
# PATH_MINIMAP2 = "/BiO/Access/yoonsung/anaconda3/envs/RagTag/bin/minimap2"


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

def run(qryfasta, reffasta, path_save, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):
    os.makedirs(os.path.dirname(path_save), exist_ok=True)

    cmd = f"{PATH_MINIMAP2} -x asm5 -o {path_save} -t {n_jobs}"
    
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

if __name__ == "__main__":
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis"
    
    curr_time = int(time.time())
    refname = "CHM13"
    path_ref = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/CHM13/chm13v2.0.fa"

    dict_hap_to_parental = {
        1:'p',
        2:'m'
    }
    
    dict_hap_to_list_region = {
        1:[
            # "chr13:1-1000000",
            # "chr13:1-20000000",
            "chrY:20000000"
        ],
        2:[
            # "chr15:98753195"
        ]
    }
    
    for additional_bin in ["Specific", "0", "a", "0a"]:
        for hap in [1, 2]:
            for region in dict_hap_to_list_region[hap]:
                region_name = region.replace(':', '_').replace('-', '_')

                # path_fasta = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/Local_Assembly/HiFi_ONT/Hap{hap}/{region_name}.with_{additional_bin}/KOREF_S1.Local_Assembly.{region_name}.ul_L50kb_QV20_Q30Perc80.triobin_{dict_hap_to_parental[hap]}_{additional_bin}.hifi_L1000_QV20.triobin_{dict_hap_to_parental[hap]}_{additional_bin}.haploid.asm.bp.p_ctg.fasta"
                # path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/Local_Assembly/HiFi_ONT/Hap{hap}/{region_name}.with_{additional_bin}/Align/Align.REF_{refname}.QRY_assembly.paf"

                path_fasta = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/Local_Assembly/ONT_only/Hap1/chrY_20000000.with_{additional_bin}/assembly.fasta"
                path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/Local_Assembly/ONT_only/Hap1/chrY_20000000.with_{additional_bin}/Align/Align.REF_{refname}.QRY_assembly.fasta"

                assert os.path.exists(path_fasta)
                run(
                    path_fasta, 
                    path_ref, 
                    path_save, 
                    f"{dir_workspace}/Minimap2_ASM_Align/tmp",
                    f"{dir_workspace}/Minimap2_ASM_Align/sh",
                    f"{dir_workspace}/Minimap2_ASM_Align/log",
                    f"Minimap2.REF_{refname}.QRY_flye_assembly.Hap{hap}.with_{additional_bin}.{region_name}.{curr_time}",
                    10,
                    [],
                    [],
                    **{}
                )
    
    

if __name__ == "__main2__":
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis"

    curr_time = int(time.time())

    modified_option = {
    }
    dict_refname_to_path = {
        "CHM13" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/CHM13/chm13v2.0.fa",
        "CN1_combine" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/CN1/CN1_combine.v1.0.1.fasta"
    }

    for refname, path_ref in dict_refname_to_path.items():
        for hap in ["1", "2"]:
            path_qry = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/Chromosome_level_scaffolds.Assign_by_CHM13/CuratedDraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.with_teloboth.hap{hap}.RagTag_with_Telomere_Assign.minimap2_latest.CHM13.pseudo_telomere_0.fasta"

            jobname = f"Align.Chromosome_Assigned.hap{hap}.Ref_{refname}"

            path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/Chromosome_level_scaffolds.Assign_by_CHM13/Minimap2_Align/CuratedDraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.with_teloboth.hap{hap}.RagTag_with_Telomere_Assign.minimap2_latest.CHM13.pseudo_telomere_0.Ref_{refname}.minimap2_2.29.paf"

            run(
                path_qry, 
                path_ref, 
                path_save, 
                f"{dir_workspace}/Minimap2_ASM_Align/tmp",
                f"{dir_workspace}/Minimap2_ASM_Align/sh",
                f"{dir_workspace}/Minimap2_ASM_Align/log",
                f"Minimap2.{jobname}.{curr_time}",
                24,
                [],
                [],
                **modified_option
            )


if __name__ == "__main2__":
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis"

    # list_path_qry_format = [
    #     (
    #         "KOREF_S1.ul_L100kb_QV15.Trio_Specific.hap{hap}.telo_both",
    #         "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/Telomere_Detection/Fasta_Split_by_Telomere/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg.telo_both.fasta",
    #         "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/Telomere_Detection/Fasta_Split_by_Telomere/Minimap2_Align/Align.KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg.telo_both"
    #     ),
    #     (
    #         "KOREF_S1.ul_L100kb_QV15.default_option.hap{hap}.telo_both",
    #         "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/Telomere_Detection/Fasta_Split_by_Telomere/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.asm.dip.hap{hap}.p_ctg.telo_both.fasta",
    #         "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/Telomere_Detection/Fasta_Split_by_Telomere/Minimap2_Align/Align.KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.asm.dip.hap{hap}.p_ctg.telo_both"
    #     )
    # ]
    # list_path_qry_format = [
    #     (
    #         "KOREF_S1.ul_L100kb_QV15.Trio_Specific.hap{hap}.ntLink.telo_both",
    #         "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling/Adjusted/Telomere_Detection/Fasta_Split_by_Telomere/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg.no_telo.telo_single_masked.fa.k32.w100.z1000.ntLink.scaffolds.adjusted_telomere.telo_both.fasta",
    #         "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling/Adjusted/Telomere_Detection/Fasta_Split_by_Telomere/Minimap2_Align/Align.KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg.no_telo.telo_single_masked.fa.k32.w100.z1000.ntLink.scaffolds.adjusted_telomere.telo_both"
    #     ),
    #     (
    #         "KOREF_S1.ul_L100kb_QV15.default_option.hap{hap}.ntLink.telo_both",
    #         "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling/Adjusted/Telomere_Detection/Fasta_Split_by_Telomere/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.asm.dip.hap{hap}.p_ctg.no_telo.telo_single_masked.fa.k32.w100.z1000.ntLink.scaffolds.adjusted_telomere.telo_both.fasta",
    #         "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling/Adjusted/Telomere_Detection/Fasta_Split_by_Telomere/Minimap2_Align/Align.KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.asm.dip.hap{hap}.p_ctg.no_telo.telo_single_masked.fa.k32.w100.z1000.ntLink.scaffolds.adjusted_telomere.telo_both"
    #     )
    # ]
    list_path_qry_format = [
        (
            "KOREF_S1.ul_L100kb_QV15.Trio_Specific.hap{hap}.ntLink.telo_single",
            "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling/Adjusted/Telomere_Detection/Fasta_Split_by_Telomere/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg.no_telo.telo_single_masked.fa.k32.w100.z1000.ntLink.scaffolds.adjusted_telomere.telo_single.fasta",
            "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling/Adjusted/Telomere_Detection/Fasta_Split_by_Telomere/Minimap2_Align.Single_Telomere.asm5/Align.KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg.no_telo.telo_single_masked.fa.k32.w100.z1000.ntLink.scaffolds.adjusted_telomere.telo_single"
        ),
        (
            "KOREF_S1.ul_L100kb_QV15.default_option.hap{hap}.ntLink.telo_single",
            "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling/Adjusted/Telomere_Detection/Fasta_Split_by_Telomere/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.asm.dip.hap{hap}.p_ctg.no_telo.telo_single_masked.fa.k32.w100.z1000.ntLink.scaffolds.adjusted_telomere.telo_single.fasta",
            "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling/Adjusted/Telomere_Detection/Fasta_Split_by_Telomere/Minimap2_Align.Single_Telomere.asm5/Align.KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.asm.dip.hap{hap}.p_ctg.no_telo.telo_single_masked.fa.k32.w100.z1000.ntLink.scaffolds.adjusted_telomere.telo_single"
        )
    ]
    
    dict_refname_to_path = {
        "CHM13" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/CHM13/chm13v2.0.fa",
        "CN1_combine" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/CN1/CN1_combine.v1.0.1.fasta"
    }
    
    curr_time = int(time.time())

    modified_option = {
    }

    for refname, path_ref in dict_refname_to_path.items():
        for hap in ["1", "2"]:
            for jobname_format, path_qry_format, path_save_prefix_format in list_path_qry_format:
                path_qry = path_qry_format.format(hap = hap)
                path_save_prefix = path_save_prefix_format.format(hap = hap)
                jobname = f"{jobname_format.format(hap = hap)}.Ref_{refname}"

                path_save = f"{path_save_prefix}.Ref_{refname}.minimap2_2.29.c_option.paf"

                run(
                    path_qry, 
                    path_ref, 
                    path_save, 
                    f"{dir_workspace}/Minimap2_ASM_Align/tmp",
                    f"{dir_workspace}/Minimap2_ASM_Align/sh",
                    f"{dir_workspace}/Minimap2_ASM_Align/log",
                    f"Minimap2.{jobname}.{curr_time}",
                    24,
                    [],
                    [],
                    **modified_option
                )
#%%