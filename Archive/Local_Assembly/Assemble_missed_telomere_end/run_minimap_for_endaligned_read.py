import os, time
import subprocess

PATH_MINIMAP2 = "/BiO/Research/Project2/Endangered_Jindo_Dog/Resources/Tools/minimap2/minimap2_v2.30_r1287/minimap2/minimap2"
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

    cmd = f"{PATH_MINIMAP2} -x map-ont -U 10,1000000 -P -c -w 10 -K 10G -t {n_jobs}"
    
    for param_name, param_value in additional_params.items():
        if len(param_name) == 1:
            cmd += f" -{param_name} {param_value}"
        else:
            cmd += f" --{param_name} {param_value}"
    cmd += f" {reffasta} {qryfasta}  | awk '{{if($7-$9<1000) print $0}}' - > {path_save}"
    
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
    dir_workspace = "/BiO/Research/Project2/Endangered_Goral/Analysis"

    curr_time = int(time.time())

    modified_option = {
    }
        
    list_contignames = [
        "h1tg000011l",
        "h1tg000012l",
        "h1tg000026l"
    ]
    
    path_ref_format = "/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/From_ONT/Reads_for_Elongate/from_contig/{contig}_elongate_1.fa"
    
    dict_qry = {
        # "Raw_HiFi":"/BiO/Research/Project2/Endangered_Goral/Results/Input_Data_Filtering/HiFi/KoGoral/Merged_Fastq/KoGoral.Merged.hifi_reads.Edge_Trim_Front_100_Tail_100.Filter_L1000_QV20.fastq.gz",
        "Raw_ONT":"/BiO/Research/Project2/Endangered_Goral/Results/Input_Data_Filtering/ONT/KoGoral/Merged_Fastq/KoGoral.Merged.ONT_reads.Adapt_Filt.Edge_Trim_Front_500_Tail_500.Filter_L50kb_QV15.fastq.gz"
    }

    for refname in list_contignames:
        for qryname, path_qry in dict_qry.items():
        # for qryname, path_qry in dict_qry.items():
            # if refname == qryname: continue
            
            path_ref = path_ref_format.format(contig = refname)
            
            jobname = f"Align.Telomere_Elongation.{qryname}.Ref_{refname}.R1.mapontPcw10U101000000"

            path_save = f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/From_ONT/RawReads_Aligned_ONT/from_round1/Minimap2_Align/Align_Minimap2.Ref_{refname}.R1.QRY_{qryname}.mapontPcw10U101000000.endaligned1000.paf"

            run(
                path_qry, 
                path_ref, 
                path_save, 
                f"{dir_workspace}/Minimap2_ASM_Align/tmp",
                f"{dir_workspace}/Minimap2_ASM_Align/sh",
                f"{dir_workspace}/Minimap2_ASM_Align/log",
                f"Minimap2.{jobname}.{curr_time}",
                256,
                [],
                [],
                **modified_option
            )

