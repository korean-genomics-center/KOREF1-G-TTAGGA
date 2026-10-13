#%%
import subprocess, os, time

PATH_SCRIPT = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Resources/Scripts/Polishing/Manual_Polishing/module_collect_lowqv_bed.py"
PATH_BEDTOOLS = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/Bedtools/bedtools"

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

def make_symlink(path_src, path_dst):
    if os.path.exists(path_dst):
        os.remove(path_dst)
    os.symlink(path_src, path_dst)

def run(path_fasta, path_qv, output_prefix, qv_thresh, extra_ranges, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list()):

    cmd = f"""
python {PATH_SCRIPT} --path_qv {path_qv} --qv_thresh {qv_thresh} --extra_ranges {extra_ranges} --out_prefix {output_prefix}.region;
{PATH_BEDTOOLS} getfasta -fi {path_fasta} -bed {output_prefix}.region.bed -fo {output_prefix}.sequence.fasta -name
{PATH_BEDTOOLS} getfasta -fi {path_fasta} -bed {output_prefix}.region.issue.bed -fo {output_prefix}.issue.sequence.fasta -name
"""
    
    os.makedirs(os.path.dirname(output_prefix), exist_ok=True)
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

    qv_thresh = 60
    extra_ranges = 250

    k = 31
    window_size = 2000
    
    n_round = 7
    afterround = 4
    
    curr_time = int(time.time())

    for hap in [1, 2]:
        path_asm = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Reference_Fasta/Manual_Polishing/After_Round{afterround}/Round{n_round-1}.Manual_Polished_Fasta/PolishedASM.KOREF_S1.Trio_Specific.UpperCase.Polish_R{n_round-1}.Man.hap{hap}.fasta"
        
        path_qv_calc = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Reference_Fasta/Manual_Polishing/After_Round{afterround}/Round{n_round-1}.Manual_Polished_Fasta/Merqury_QV_per_window/PolishedASM.KOREF_S1.Trio_Specific.UpperCase.Polish_R{n_round-1}.Man.hap{hap}.QV_per_window.k{k}.chunked.{window_size}.qv"

        output_prefix = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Manual_Polishing/Round{n_round}_afterR{afterround}/Collect_LowQV_Region/LowQVRegion.KOREF_S1.Trio_Specific.UpperCase.Polish_R{n_round-1}.Man.hap{hap}.QV_per_window.k{k}.chunked.{window_size}.qv_less{qv_thresh}.extra_range{extra_ranges}"

        assert os.path.exists(path_asm)
        assert os.path.exists(path_qv_calc)

        run(
            path_asm,
            path_qv_calc,
            output_prefix,
            qv_thresh,
            extra_ranges,
            f"{dir_workspace}/Polishing/Manual_Polishing/Collect_LowQV/tmp",
            f"{dir_workspace}/Polishing/Manual_Polishing/Collect_LowQV/sh",
            f"{dir_workspace}/Polishing/Manual_Polishing/Collect_LowQV/log",
            f"Collect_LowQV.CuratedGapclosedASM.KOREF_S1.Polish_R{n_round-1}.Hap{hap}.k{k}.{window_size}.lowqv_{qv_thresh}.extra_range{extra_ranges}.{curr_time}",
            4,
            [],
            []
        )
#%%