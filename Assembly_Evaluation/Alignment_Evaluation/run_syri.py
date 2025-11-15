import os, time
import subprocess

PATH_SYRI = "/BiO/Access/yoonsung/anaconda3/envs/syri/bin/syri"
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

def run(path_paf, reffasta, qryfasta, path_output_prefix, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):
    dir_save = os.path.dirname(path_output_prefix)
    os.makedirs(dir_save, exist_ok=True)
    
    make_symlink(path_paf, dir_save)
    make_symlink(reffasta, dir_save)
    make_symlink(qryfasta, dir_save)

    cmd = f"cd {dir_save};{PATH_SYRI} -c {os.path.basename(path_paf)} -r {os.path.basename(reffasta)} -q {os.path.basename(qryfasta)} -F P --dir {dir_save} --prefix {os.path.basename(path_output_prefix)} --nc {n_jobs} --no-chrmatch"
    
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
    
    dict_refname_to_path = {
        "CHM13" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/CHM13/chm13v2.0.fa",
        "CN1_Pat" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/CN1/CN1_pat.v1.0.1.fasta",
        "CN1_Mat" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/CN1/CN1_mat.v1.0.1.fasta",
        "YAO_Pat" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/T2T_YAO/T2T_Yao_pat_v1.1.fasta",
        "YAO_Mat" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/T2T_YAO/T2T_Yao_mat_v1.1.fasta",
        "I002C_Pat" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/I002C/I002C_Paternal_v0.7.fasta",
        "I002C_Mat" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/I002C/I002C_Maternal_v0.7.fasta",
        "KOREF_Hap1":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Final_Assembly/FinalASM.KOREF_S1.Trio_Specific.hap1.Reordered.fasta",
        "KOREF_Hap2":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Final_Assembly/FinalASM.KOREF_S1.Trio_Specific.hap2.Reordered.fasta"
    }

    for (refname, qryname) in [("CHM13", "KOREF_Hap2")]:
        
        reffasta = dict_refname_to_path[refname]
        qryfasta = dict_refname_to_path[qryname]
        path_paf = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Alignment_Evaluation/Align_Minimap2/Minimap2_Align.Ref_{refname}.QRY_{qryname}.asm5c.paf"

        jobname = f"REF_{refname}.QRY_{qryname}"

        path_output_prefix = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Alignment_Evaluation/SyRI/{refname}.{qryname}/syri.ref_{refname}.qry_{qryname}"

        assert os.path.exists(reffasta)
        assert os.path.exists(qryfasta)
        assert os.path.exists(path_paf)
        
        run(
            path_paf,
            reffasta,
            qryfasta,
            path_output_prefix, 
            f"{dir_workspace}/Alignment_Evaluation/tmp",
            f"{dir_workspace}/Alignment_Evaluation/sh",
            f"{dir_workspace}/Alignment_Evaluation/log",
            f"Run_SyRI.{jobname}.{curr_time}",
            16,
            [],
            [],
            **modified_option
        )
#%%