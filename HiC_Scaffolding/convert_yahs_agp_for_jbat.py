#%%
import subprocess, os, time

PATH_YAHS_JUICER = "/BiO/Research/Project2/Endangered_Goral/Resources/Tools/YaHS/YaHS-v1.2.2/yahs/juicer"
PATH_JUICER_TOOLS = "/BiO/Research/Project2/Endangered_Goral/Resources/Tools/Juicer/Juicer_v1.6_master/juicer_1.6_CPU/juicer_tools_1.22.01.jar"


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

def run(path_bin, path_agp, path_fai, path_output_prefix, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):

    cmd = f"""
{PATH_YAHS_JUICER} pre -a -o {path_output_prefix}.pre {path_bin} {path_agp} {path_fai} > {path_output_prefix}.pre.log 2>&1;
(java -jar -Xms128G -Xmx128G {PATH_JUICER_TOOLS} pre {path_output_prefix}.pre.txt {path_output_prefix}.hic.part <(cat {path_output_prefix}.pre.log  | grep PRE_C_SIZE | awk '{{print $2" "$3}}')) && (mv {path_output_prefix}.hic.part {path_output_prefix}.hic)
"""
    
    os.makedirs(os.path.dirname(path_output_prefix), exist_ok=True)
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
    
    for hap in [1]:
        path_bin = "/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/HiC_Scaffolding_with_additional_contigs/YaHS_Run.with_additional_teloseq.with_raw_unclustered_read/YaHS_Run.KoGoral_DraftASM.Hap1.with_additional_teloseq.with_raw_unclustered_read.bin"
        path_agp = "/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/HiC_Scaffolding_with_additional_contigs/YaHS_Run.with_additional_teloseq.with_raw_unclustered_read/YaHS_Run.KoGoral_DraftASM.Hap1.with_additional_teloseq.with_raw_unclustered_read_scaffolds_final.agp"
        path_fai = "/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/HiC_Scaffolding_with_additional_contigs/HiC_Alignment.4Enz.with_raw_unclustered_read/genome/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.asm.hic.hap1.p_ctg.with_additional_teloseq.with_unclustered_rawONT.fasta.fai"
        
        path_output_prefix = f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/HiC_Scaffolding_with_additional_contigs/YaHS_Run.with_additional_teloseq.with_raw_unclustered_read/JBAT/JBAT_out.KoGoral_DraftASM.Hap1.with_additional_teloseq.with_raw_unclustered_read"
        
        curr_time = int(time.time())
        run(
            path_bin, path_agp, path_fai, path_output_prefix,
            f"{dir_workspace}/Convert_YaHS_to_HiC/tmp",
            f"{dir_workspace}/Convert_YaHS_to_HiC/sh",
            f"{dir_workspace}/Convert_YaHS_to_HiC/log",
            f"Convert_YaHS_to_HiC.KoGoral.DraftASM.Hap{hap}.with_additional_teloseq.with_raw_unclustered_read.{curr_time}",
            8,
            [],
            []
        )
#%%