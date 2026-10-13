import os, time
import subprocess
import pandas as pd

PATH_PLOTSR = "/BiO/Access/yoonsung/anaconda3/envs/syri/bin/plotsr"
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

def run(list_syri_out, list_asmnames, list_asmpaths, list_chr_order, path_output_prefix, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):
    dir_save = os.path.dirname(path_output_prefix)
    os.makedirs(dir_save, exist_ok=True)
    
    table_asm = pd.DataFrame({
        "#file":list_asmpaths,
        "name":list_asmnames
    })
    
    # path_save_chrorder = f"{path_output_prefix}.chrorder.txt"
    path_save_genomes = f"{path_output_prefix}.genomes.txt"
    
    table_asm.to_csv(path_save_genomes, sep = '\t', index = False)
    # with open(path_save_chrorder, 'w') as fw:
    #     for chrname in list_chr_order:
    #         fw.write(f"{chrname}\n")
    
    syri_inputs = ' '.join(list(map(lambda path_syri: f"--sr {path_syri}", list_syri_out)))
    chrs_select = ' '.join(list(map(lambda chrname: f"--chr {chrname}", list_chr_order)))

    cmd = f"{PATH_PLOTSR} {syri_inputs} --genomes {path_save_genomes} {chrs_select} -o {path_output_prefix}.pdf "
    
    for key, val in additional_params.items():
        cmd += f"--{key} {val} "
    
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

    modified_option = {
    }
    
    list_syri_out = [
        "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Alignment_Evaluation/SyRI/KOREF_Hap1.CHM13/syri.out",
        "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Alignment_Evaluation/SyRI/CHM13.KOREF_Hap2/syri.ref_CHM13.qry_KOREF_Hap2syri.out"
    ]
    
    list_asmnames = [
        "KOREF_Hap1", "CHM13", "KOREF_Hap2"
    ]
    
    list_asmpaths = [
        "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Final_Assembly/FinalASM.KOREF_S1.Trio_Specific.hap1.Reordered.fasta",
        "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/CHM13/chm13v2.0.fa",
        "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Final_Assembly/FinalASM.KOREF_S1.Trio_Specific.hap2.Reordered.fasta"
    ]
    
    list_chr_order = list(map(lambda val: f"chr{val}", range(1, 23))) # + ["chrX", "chrY"]
    
    path_output_prefix = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Alignment_Evaluation/PlotSR/{'_'.join(list_asmnames)}/plotsr.{'_'.join(list_asmnames)}"
    
    jobname = f"Syri_{'_'.join(list_asmnames)}"
    
    run(
        list_syri_out,
        list_asmnames,
        list_asmpaths,
        list_chr_order,
        path_output_prefix, 
        f"{dir_workspace}/Alignment_Evaluation/tmp",
        f"{dir_workspace}/Alignment_Evaluation/sh",
        f"{dir_workspace}/Alignment_Evaluation/log",
        f"Run_plotSR.{jobname}.{curr_time}",
        16,
        [],
        [],
        **modified_option
    )
#%%