#%%
import subprocess, os, time

PATH_RAGTAG_SCAFFOLD = "/BiO/Access/yoonsung/anaconda3/envs/RagTag/bin/ragtag_scaffold_asm20.py"
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

def write_shell_script_and_run(path_shell, shell_context, env = None):
    with open(path_shell, 'w') as fw:
        fw.write(shell_context)
    
    subprocess.run(f"qsub {path_shell}", shell = True, env = env)

def run(reference_fasta, scaffolds_fasta, path_dir_save, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):
    cmd = f"""{PATH_RAGTAG_SCAFFOLD} -r -o {path_dir_save} -w -t {n_jobs} --debug --aligner {PATH_MINIMAP2} {reference_fasta} {scaffolds_fasta}"""
    
    for param_name, param_value in additional_params.items():
        if len(param_name) == 1:
            cmd += f" {param_name}={param_value}"
        else:
            cmd += f" {param_name}={param_value}"
    
    # merqury_environ = os.environ.copy()
    # merqury_environ["PATH"] = f"{DIR_MERYL}:{DIR_SAMTOOLS}:{DIR_BEDTOOLS}:{merqury_environ['PATH']}"
    # merqury_environ["MERQURY"] = DIR_MERQURY
    # merqury_environ["OMP_NUM_THREADS"] = f"{n_jobs}"
    
    os.makedirs(path_dir_save, exist_ok=True)
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
        env = None
    )

if __name__ == "__main__":
    
    dict_hap_to_side = {
        "1" : "p",
        "2" : "m"
    }
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis"

    dict_refname_to_path = {
        "CHM13" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/CHM13/chm13v2.0.fa",
        "CN1_combine" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/CN1/CN1_combine.v1.0.1.fasta"
    }
    for input_fasta_type in ["no_telo.telo_single", "telo_single"]:
        # dict_refname_to_path = {
        #     "CHM13" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling.RagTag.default_option.CHM13/CHM13.without_assigned_chromosomes_hap{hap}.fasta",
        #     "CN1_combine" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling.RagTag.default_option.CN1_combine/CN1_combine.without_assigned_chromosomes_hap{hap}.fasta"
        # }
        for hap in ["1", "2"]:
            for refname in dict_refname_to_path.keys():
                reference_fasta = dict_refname_to_path[refname]
                scaffolds_fasta = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling/Adjusted/Telomere_Detection/Fasta_Split_by_Telomere/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.asm.dip.hap{hap}.p_ctg.no_telo.telo_single_masked.fa.k32.w100.z1000.ntLink.scaffolds.adjusted_telomere.{input_fasta_type}.fasta"
                
                path_dir_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling.RagTag.default_option.{refname}/DraftASM.KOREF_S1.ul_L100kb_QV15.default_option.ntLink.adjusted_telomere.{input_fasta_type}.hap{hap}.RagTag.minimap_asm20.minimap2_latest.{refname}"
                
                curr_time = int(time.time())
                run(
                    reference_fasta, 
                    scaffolds_fasta, 
                    path_dir_save,
                    f"{dir_workspace}/RagTag_Scaffolding/tmp",
                    f"{dir_workspace}/RagTag_Scaffolding/sh",
                    f"{dir_workspace}/RagTag_Scaffolding/log",
                    f"Scaffold_RagTag.KOREF_ASM.ul_L100kb_QV15.default_option.ntLink.adjusted_telomere.{input_fasta_type}.hap{hap}.{refname}.{curr_time}",
                    16,
                    [],
                    []
                )

    for input_fasta_type in ["no_telo.telo_single", "telo_single"]:
        # dict_refname_to_path = {
        #     "CHM13" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling.RagTag.default_option.CHM13/CHM13.without_assigned_chromosomes_hap{hap}.fasta",
        #     "CN1_combine" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling.RagTag.default_option.CN1_combine/CN1_combine.without_assigned_chromosomes_hap{hap}.fasta"
        # }
        for hap in ["1", "2"]:
            for refname in dict_refname_to_path.keys():
                reference_fasta = dict_refname_to_path[refname]
                scaffolds_fasta = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling/Adjusted/Telomere_Detection/Fasta_Split_by_Telomere/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg.no_telo.telo_single_masked.fa.k32.w100.z1000.ntLink.scaffolds.adjusted_telomere.{input_fasta_type}.fasta"
                
                path_dir_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling.RagTag.default_option.{refname}/DraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.{input_fasta_type}.hap{hap}.RagTag.minimap_asm20.minimap2_latest.{refname}"
                
                curr_time = int(time.time())
                run(
                    reference_fasta, 
                    scaffolds_fasta, 
                    path_dir_save,
                    f"{dir_workspace}/RagTag_Scaffolding/tmp",
                    f"{dir_workspace}/RagTag_Scaffolding/sh",
                    f"{dir_workspace}/RagTag_Scaffolding/log",
                    f"Scaffold_RagTag.KOREF_ASM.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.{input_fasta_type}.hap{hap}.{refname}.{curr_time}",
                    16,
                    [],
                    []
                )
#%%