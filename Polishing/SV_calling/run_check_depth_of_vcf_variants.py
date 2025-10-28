import os, time
import subprocess

PATH_SCRIPT = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Resources/Scripts/Polishing/SV_calling/check_depth_of_vcf_variants.py"


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

def run(path_vcf, path_mosdepth, path_output, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):

    cmd = f"python {PATH_SCRIPT} --vcf {path_vcf} --mosdepth {path_mosdepth} --out {path_output} --n_jobs {n_jobs}"
    
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

    dict_tgs_to_name = {
        "HiFi":"PacBio_HiFi_ccs_converted.Filter_L1000_QV20",
        "ONT":"ONT_UL.Filter_L10kb_QV15"
    }
    dict_hap_to_parental = {
        1:"p",
        2:"m"
    }
    
    curr_time = int(time.time())
    
    for hap in [1, 2]:
        path_vcf = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round8_SV_Call/SV_Call_with_Unclassified/SV_Calling/Merged_SV/SV_Call.KOREF_PolishedASM_Round8_SV.hap{hap}.hifi_ont.sametrio.filtered.iris.merged.split.homoalt_include.hetalt_only.filt_difftrio_homoref.filt_no_asmsv.vaf_0.5.minsupp_5.call_merged.indel_only.vcf"
        
        for tgs_type, tgs_name in dict_tgs_to_name.items():
            path_mosdepth = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round8_SV_Call/SV_Call_without_Unclassified/Alignment_Depth/Mosdepth.Primary.Ref.KOREF_PolishedASM.Polish_R7.Man.hap{hap}.QRY.{tgs_name}.triobin_split.triobin_{dict_hap_to_parental[hap]}.per-base.bed.gz"
            
            path_output = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round8_SV_Call/SV_Call_Merged/Aligned_Depth_on_TrioSpecific/Aligned_Depth.hap{hap}.{tgs_name}.Trio_Specific_only.triobin_{dict_hap_to_parental[hap]}.of_SVs_with_Unclassified.tsv"
            
            assert os.path.exists(path_vcf), path_vcf
            assert os.path.exists(path_mosdepth), path_mosdepth

            run(
                path_vcf, 
                path_mosdepth, 
                path_output,
                f"{dir_workspace}/Polishing/Polishing_SV/Check_Depth_on_SV/tmp",
                f"{dir_workspace}/Polishing/Polishing_SV/Check_Depth_on_SV/sh",
                f"{dir_workspace}/Polishing/Polishing_SV/Check_Depth_on_SV/log",
                jobname = f"Check_Depth_on_SV.hap{hap}.SV_with_Unclassified.On_{tgs_type}_triobin_{dict_hap_to_parental[hap]}.{curr_time}",
                n_jobs = 30
            )