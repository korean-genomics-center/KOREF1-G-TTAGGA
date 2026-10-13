import pandas as pd
import os, time
import subprocess

PATH_KENTTOOOLS = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/CAT/Dependencies/kent"
PATH_HALTOOLS = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/CAT/Dependencies/cactus/cactus-bin-v3.1.0/bin"

PATH_AUGUSTUS_CONFIG = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/CAT/CAT_v2.2.1/Comparative-Annotation-Toolkit/augustus_cfgs"


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

def run(path_hal, list_path_shortread_bam, list_path_isoseq_bam, path_ref_gff3, refname, asmname, dir_out, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list()):    
    run_cgp = False
    run_pb = False

    dir_out_work = os.path.join(dir_out, "work_dir")
    os.makedirs(dir_out, exist_ok=True)
    os.makedirs(dir_out_work, exist_ok=True)

    path_save_config = os.path.join(dir_out, "runcat.config")
    with open(path_save_config, 'w') as fw:
        fw.write(f"[ANNOTATION]\n{refname} = {path_ref_gff3}\n\n")

        if len(list_path_shortread_bam) > 0:
            fw.write(f"[BAM]\n{asmname} = {', '.join(list_path_shortread_bam)}\n\n")
            run_cgp = True
            
        if len(list_path_isoseq_bam) > 0:
            fw.write(f"[ISO_SEQ_BAM]\n{asmname} = {', '.join(list_path_isoseq_bam)}\n\n")
            run_pb = True

    cmd = f"""
luigi --module cat RunCat --hal={path_hal} --ref-genome={refname} --target-genomes='("{asmname}",)' --workers={n_jobs} --maxCores {n_jobs} --config={path_save_config} --work-dir {dir_out_work} --out-dir {dir_out} --local-scheduler --augustus {'--augustus-cgp' if run_cgp else ''} {'--augustus-pb' if run_pb else ''} --assembly-hub --binary-mode local --logLevel INFO --tm-cfg {PATH_AUGUSTUS_CONFIG}/extrinsic.ETM1.cfg --tmr-cfg {PATH_AUGUSTUS_CONFIG}/extrinsic.ETM2.cfg --augustus-cgp-cfg-template {PATH_AUGUSTUS_CONFIG}/cgp_extrinsic_template.cfg --pb-cfg {PATH_AUGUSTUS_CONFIG}/extrinsic.M.RM.PB.E.W.cfg
"""
    
    runcat_env = os.environ
    runcat_env["PATH"] = f"""{runcat_env["PATH"]}:{PATH_KENTTOOOLS}:{PATH_HALTOOLS}"""
        
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


    for hap in [1, 2]:
        path_hal = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Cactus_Alignment/Align.Hap{hap}/grch38_koref1hap{hap}.full.hal"

        list_path_shortread_bam = [
            f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/ShortRead_Alignment/STAR_Alignment/Hap{hap}/STAR_Align.Ref_KOREF1_FinalASM_Hap{hap}.QRY_Blood_Trim_F15_T5.twopassModeAligned.sortedByCoord.out.sorted.bam",
            f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/ShortRead_Alignment/STAR_Alignment/Hap{hap}/STAR_Align.Ref_KOREF1_FinalASM_Hap{hap}.QRY_CellLine_Trim_F15_T5.twopassModeAligned.sortedByCoord.out.sorted.bam"
        ]
        list_path_isoseq_bam = [
            f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/IsoSeq_Preprocessing/Isoseq_align/Minimap2_align.Ref_KOREF1_Hap{hap}.Qry_Isoseq_all.splicehqMDeqxNosec.sorted.bam"
        ]

        path_ref_gff3 = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Input_GFF_Processing/gencode.v49.annotation.NoPARY.NochrM.gff3"

        refname = "GRCh38"
        asmname = f"KOREF1_ASM_Hap{hap}"

        dir_out = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Run_CAT/Hap{hap}.Blood_Cellline.CGP_PB"

        run(
            path_hal, list_path_shortread_bam, list_path_isoseq_bam, path_ref_gff3, refname, asmname, dir_out, 
            f"{dir_workspace}/CAT_Annotation/Run_Cat/tmp", 
            f"{dir_workspace}/CAT_Annotation/Run_Cat/sh", 
            f"{dir_workspace}/CAT_Annotation/Run_Cat/log", 
            f"Run_Cat.Ref.KOREF1_FinalASM.hap{hap}.Blood_Cellline.CGP_PB.{curr_time}",
            96,
            jobs_wait=[]
        )