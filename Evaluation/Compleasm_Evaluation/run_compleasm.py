import os, time
import subprocess

PATH_COMPLEASM = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/Compleasm/Compleasm_v0.2.7/compleasm_kit/compleasm.py"
PATH_HMMSEARCH = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/Compleasm/Compleasm_v0.2.7/compleasm_kit/hmmsearch"
PATH_MINIPROT = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/Compleasm/Compleasm_v0.2.7/compleasm_kit/miniprot"


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

def run(path_assembly, dir_output, lineage, path_library, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list()):    

    cmd = f"""python {PATH_COMPLEASM} run -a {path_assembly} -o {dir_output} -t {n_jobs} -l {lineage} -L {path_library} --miniprot_execute_path {PATH_MINIPROT} --hmmsearch_execute_path {PATH_HMMSEARCH}"""
        
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

    path_library = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/Compleasm/Compleasm_v0.2.7/compleasm_kit/mb_downloads"

    dict_asm_to_files = {
        "CHM13_v2":{
            "Haploid":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/CHM13/chm13v2.0.fa"
        },
        "T2T_YAO":{
            "Paternal":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/T2T_YAO/T2T_Yao_pat_v1.1.fasta",
            "Maternal":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/T2T_YAO/T2T_Yao_mat_v1.1.fasta",
            "Haploid":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/T2T_YAO/T2T_Yao_hp_v1.1.fasta"
        },
        "CN1":{
            "Paternal":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/CN1/CN1_pat.v1.0.1.fasta",
            "Maternal":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/CN1/CN1_mat.v1.0.1.fasta",
            "Haploid":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/CN1/CN1_combine.v1.0.1.fasta"
        },
        "I002C":{
            "Paternal" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/I002C/I002C_Paternal_v0.7.fasta",
            "Maternal" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/I002C/I002C_Maternal_v0.7.fasta"
        }
    }
    
    for asmname in dict_asm_to_files.keys():
        for partname, path_assembly in dict_asm_to_files[asmname].items():
            dir_output = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Compleasm_Evaluation/{asmname}/{partname}"

            assert os.path.exists(path_assembly)
            assert os.path.exists(path_library)
            
            lineage = "primates"

            
            curr_time = int(time.time())
            run(
                path_assembly,
                dir_output,
                lineage,
                path_library,
                f"{dir_workspace}/Compleasm_Evaluation/Run_Compleasm/tmp",
                f"{dir_workspace}/Compleasm_Evaluation/Run_Compleasm/sh",
                f"{dir_workspace}/Compleasm_Evaluation/Run_Compleasm/log",
                f"Run_Compleasm.{asmname}.{partname}.Autosome.{curr_time}",
                8,
                hostname = [],
                jobs_wait = []
            )