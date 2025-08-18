import pandas as pd
import os, time
import subprocess

PATH_WINNNOWMAP = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/winnowmap/winnowmap.v2.03_fork/Winnowmap/bin/winnowmap"
PATH_SAMTOOLS = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/Samtools/samtools-1.17/samtools-1.17/bin/samtools"


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

def run(path_fasta_qry, path_fasta_ref, path_ref_merylcount, qry_type, path_save, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list()):    

    cmd = f"""{PATH_WINNNOWMAP} --MD -W {path_ref_merylcount} -t {n_jobs} -ax {qry_type} -K 40G {path_fasta_ref} {path_fasta_qry} | {PATH_SAMTOOLS} view -bS -@ {min(n_jobs, 16)} -o {path_save}"""
        
    os.makedirs(os.path.dirname(path_save), exist_ok=True)
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

    path_pacbio_file_format = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Read_Alignment/Parental_Split_Data/PacBio_HiFi_ccs_converted.filt.Filter_L1000_QV20.merged.triobin_split.triobin_{parental}.fq.gz"
    path_ont_file_format = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Read_Alignment/Parental_Split_Data/ONT_UL.Adapt_Filt.Edge_Trim_Front_1500_Tail_1500.Filter_L10kb_QV15.merged.triobin_split.triobin_{parental}.fq.gz"
        
    dir_save = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Read_Alignment/Align.Winnowmap2/Align_Polished_ASM.From_R2"

    list_qry_fasta = [
        ("map-ont", path_ont_file_format, "ONT_UL.Filter_L10kb_QV15"),
        ("map-pb", path_pacbio_file_format, "PacBio_HiFi_ccs_converted.Filter_L1000_QV20")
    ]

    for qry_type, path_qry_format, jobname_suffix in list_qry_fasta:
        for hap in [1, 2]:
            
            path_fasta_ref = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Reference_Fasta/Round2.SNV_Polished_Fasta/PolishedASM.KOREF_S1.Trio_Specific.UpperCase.Polish_R2.SNV.hap{hap}.fasta"
            path_ref_merylcount = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Read_Alignment/MerylDB/KOREF_PolishedASM.R2_SNV.Hap{hap}.kmer15.repetitive_k15.txt"

            for parental in ['m', 'p', 'a', '0']:
                    path_fasta_qry = path_qry_format.format(parental = parental)

                    path_save = os.path.join(dir_save, f"Winnowmap_Align.Ref.KOREF_PolishedASM.R2_SNV.hap{hap}.QRY.{jobname_suffix}.triobin_split.triobin_{parental}.bam")

                    assert os.path.exists(path_fasta_qry)
                    assert os.path.exists(path_fasta_ref)
                    assert os.path.exists(path_ref_merylcount)

                    curr_time = int(time.time())
                    run(
                        path_fasta_qry, 
                        path_fasta_ref, 
                        path_ref_merylcount, 
                        qry_type, 
                        path_save,
                        f"{dir_workspace}/Polishing/Winnowmap_Align/tmp",
                        f"{dir_workspace}/Polishing/Winnowmap_Align/sh",
                        f"{dir_workspace}/Polishing/Winnowmap_Align/log",
                        f"Winnowmap_Align.{jobname_suffix}.REF.KOREF_PolishedASM.R2_SNV.hap{hap}.triobin_{parental}.{curr_time}",
                        64,
                        hostname = [],
                        jobs_wait = []
                    )