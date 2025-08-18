import pandas as pd
import os, time
import subprocess

PATH_MERFIN = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/Merfin/Merfin_v1.1/merfin/build/bin/merfin"


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

def run(path_reference, path_readmer, path_lookup, path_vcf, peak, path_output_prefix, dir_tmp, dir_shell, dir_log, jobname, n_jobs, memory = 100, hostname = list(), jobs_wait = list()):   
    dir_out = os.path.dirname(path_output_prefix)

    cmd = f"""cd {dir_out};{PATH_MERFIN} -polish -sequence {path_reference} -readmers {path_readmer} -peak {peak} -prob {path_lookup} -vcf {path_vcf} -output {path_output_prefix} -threads {n_jobs} -memory {memory} -debug"""
    
    os.makedirs(dir_out, exist_ok=True)
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

    path_readmer = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Read_Alignment/MerylDB/GenomeScope2/PacBio_HiFi_ccs_converted.filt.Filter_L1000_QV20.k21.readmer.gt1.meryl"
    # path_lookup = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Read_Alignment/MerylDB/GenomeScope2/PacBio_HiFi_ccs_converted.filt.Filter_L1000_QV20.k21.readmer.genomescope_out/lookup_table.txt"
    # peak = 88.29
    path_lookup = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Read_Alignment/MerylDB/GenomeScope2/PacBio_HiFi_ccs_converted.filt.Filter_L1000_QV20.k21.l40.readmer.genomescope_out/lookup_table.txt"
    peak = 44.2

    kmer_size = 21
    
    dict_hap_to_parental = {
        1:'p',
        2:'m'
    }

    
    curr_time = int(time.time())

    for hap in [1, 2]:
        path_reference = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Reference_Fasta/Gapclosed_Fasta/CuratedGapclosedASM.KOREF_S1.Trio_Specific.UpperCase.hap{hap}.fasta"
        path_vcf = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Align_Gapclosed_ASM/SV_Call.Inc_Amb/SV_Calling/Merged_SV/SV_Call.KOREF_GapclosedASM_SV_Polishing.Inc_Amb.hap{hap}.hifi_ont.sametrio.filtered.iris.merged.split.homoalt_include.hetalt_only.filt_difftrio_homoref.filt_no_asmsv.vaf_0.5.minsupp_5.call_merged.iris.vcf"
        path_output_prefix = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Align_Gapclosed_ASM/SV_Call.Inc_Amb/Merfin_Filtered/SV_Call.KOREF_GapclosedASM_SV_Polishing.Inc_Amb.hap{hap}.hifi_ont.sametrio.filtered.iris.merged.split.homoalt_include.hetalt_only.filt_difftrio_homoref.filt_no_asmsv.vaf_0.5.minsupp_5.call_merged.iris"

        name = f"SVcall.Hap{hap}.From_Sniffles2.hifi_ont.call_merged.iris.k{kmer_size}"

        assert os.path.exists(path_reference)
        assert os.path.exists(path_vcf)
        
        curr_time = int(time.time())
        run(
            path_reference,
            path_readmer,
            path_lookup,
            path_vcf,
            peak,
            path_output_prefix,
            f"{dir_workspace}/Polishing/Polishing_SV/Merfin_Filtering/tmp",
            f"{dir_workspace}/Polishing/Polishing_SV/Merfin_Filtering/sh",
            f"{dir_workspace}/Polishing/Polishing_SV/Merfin_Filtering/log",
            f"Merfin_Filtering.{name}.{curr_time}",
            32,
            memory = 1000,
            hostname = [],
            jobs_wait = []
        )