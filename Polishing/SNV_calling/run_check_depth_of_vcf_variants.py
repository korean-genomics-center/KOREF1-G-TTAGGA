import os, time
import subprocess

PATH_SCRIPT = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Resources/Scripts/Polishing/SNV_calling/check_depth_of_vcf_variants.py"


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
    

if __name__ == "__main2__":
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
        for tgs_type, tgs_name in dict_tgs_to_name.items():
            path_vcf = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round8_SNV_Call/SNV_Filtering/Triobinned_for_Evaluation/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.Polish_R8.SV.hap{hap}.{tgs_type}.triobin_{dict_hap_to_parental[hap]}0a.PASS.NonHomRef.Multiallelic_Norm.VAF05.NonHomRef.vcf.gz"
            
            path_mosdepth = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round8_SNV_Call/SNV_Filtering/Alignment_Depth/Mosdepth.Primary.Ref.KOREF_PolishedASM.Polish_R8.SV.hap{hap}.QRY.{tgs_name}.triobin_split.triobin_{dict_hap_to_parental[hap]}.per-base.bed.gz"
            
            path_output = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round8_SNV_Call/SNV_Filtering/Triobinned_for_Evaluation/Aligned_Depth.hap{hap}.{tgs_type}.triobin_{dict_hap_to_parental[hap]}0a.PASS.NonHomRef.Multiallelic_Norm.VAF05.NonHomRef.on_triobin_{dict_hap_to_parental[hap]}.tsv"
            
            assert os.path.exists(path_vcf), path_vcf
            assert os.path.exists(path_mosdepth), path_mosdepth

            run(
                path_vcf, 
                path_mosdepth, 
                path_output,
                f"{dir_workspace}/Polishing/Polishing_SNV/Check_Depth_on_SNV/tmp",
                f"{dir_workspace}/Polishing/Polishing_SNV/Check_Depth_on_SNV/sh",
                f"{dir_workspace}/Polishing/Polishing_SNV/Check_Depth_on_SNV/log",
                jobname = f"Check_Depth_on_SNV.hap{hap}.{tgs_type}.triobin_{dict_hap_to_parental[hap]}0a.from_triobin_{dict_hap_to_parental[hap]}.{curr_time}",
                n_jobs = 30
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
        path_vcf = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round8_SNV_Call/SNV_Filtering/Triobinned_ONT_HomoAlt/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.Polish_R8.SV.hap{hap}.ONT.triobin_{dict_hap_to_parental[hap]}0a.PASS.NonHomRef.biallelic.VAF08.DP10.Norm.vcf.gz"
        
        path_mosdepth = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round8_SNV_Call/SNV_Filtering/Alignment_Depth/Mosdepth.Primary.Ref.KOREF_PolishedASM.Polish_R8.SV.hap{hap}.QRY.ONT_UL.Filter_L10kb_QV15.triobin_split.triobin_{dict_hap_to_parental[hap]}.per-base.bed.gz"
            
        path_output = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round8_SNV_Call/SNV_Filtering/Triobinned_ONT_HomoAlt/Aligned_Depth.hap{hap}.ONT.triobin_{dict_hap_to_parental[hap]}0a.PASS.NonHomRef.biallelic.VAF08.DP10.Norm.on_triobin_{dict_hap_to_parental[hap]}.tsv"
        
        assert os.path.exists(path_vcf), path_vcf
        assert os.path.exists(path_mosdepth), path_mosdepth

        run(
            path_vcf, 
            path_mosdepth, 
            path_output,
            f"{dir_workspace}/Polishing/Polishing_SNV/Check_Depth_on_SNV/tmp",
            f"{dir_workspace}/Polishing/Polishing_SNV/Check_Depth_on_SNV/sh",
            f"{dir_workspace}/Polishing/Polishing_SNV/Check_Depth_on_SNV/log",
            jobname = f"Check_Depth_on_SNV.hap{hap}.Triobinned_ONT_HomoAlt.triobin_{dict_hap_to_parental[hap]}0a.from_triobin_{dict_hap_to_parental[hap]}.{curr_time}",
            n_jobs = 30
        )