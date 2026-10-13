#%%
import subprocess, os, time

PATH_GFFREAD = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/gffread-0.12.7.Linux_x86_64/gffread"

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

def run(path_genome, path_gff, path_output_prefix, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):
    dir_output = os.path.dirname(path_output_prefix)

    cmd = f"""
{PATH_GFFREAD} -g {path_genome} -w {path_output_prefix}.transcript.fa {path_gff};
{PATH_GFFREAD} -g {path_genome} -x {path_output_prefix}.transcript_CDS.fa {path_gff};
{PATH_GFFREAD} -g {path_genome} -y {path_output_prefix}.protein.fa -S {path_gff}
{PATH_GFFREAD} -g {path_genome} -y {path_output_prefix}.protein.skip_internalstop.fa -J {path_gff}
"""
    
    os.makedirs(dir_output, exist_ok=True)
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

    dict_gff_path_format = {
        "whole_gff":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Run_CAT/Hap{hap}.Blood_Cellline.CGP_PB/consensus_gene_set/KOREF1_ASM_Hap{hap}.gff3",
        "unknown_gff":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/Process_Unknown/Hap{hap}.Blood_Cellline.CGP_PB/KOREF1_ASM_Hap{hap}.unknown_full.gff3",
        "known_gff":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/Process_Unknown/Hap{hap}.Blood_Cellline.CGP_PB/KOREF1_ASM_Hap{hap}.known_full.gff3"
    }

    path_output_prefix_format = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/CDS_Fasta/Hap{hap}.Blood_Cellline.CGP_PB/KOREF1_ASM_Hap{hap}.{gff_type}"

    for hap in [1, 2]:
        path_genome = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Final_Assembly/FinalASM.KOREF_S1.Trio_Specific.hap{hap}.Reordered.fasta"

        for gff_type, path_gff_format in dict_gff_path_format.items():
            path_gff = path_gff_format.format(hap = hap)
            path_output_prefix = path_output_prefix_format.format(hap = hap, gff_type = gff_type)
            
            run(
                path_genome, path_gff, path_output_prefix,
                f"{dir_workspace}/GFF_Process/tmp",
                f"{dir_workspace}/GFF_Process/sh",
                f"{dir_workspace}/GFF_Process/log",
                f"Generate_CDS_fasta.hap{hap}.{gff_type}.{curr_time}",
                1,
                [],
                []
            )