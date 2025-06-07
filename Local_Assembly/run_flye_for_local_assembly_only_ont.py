import pandas as pd
import os, time
import subprocess

PATH_FLYE = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/Flye/Flye_v2.9.6/Flye/bin/flye"


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

def run(path_ont_hq, dir_out, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list()):    

    cmd = f"""{PATH_FLYE} --nano-raw {path_ont_hq} -o {dir_out} -t {n_jobs}"""
        
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
    
    dict_hap_to_parental = {
        1:'p',
        2:'m'
    }
    
    dict_name_to_triobins = {
        1:{
            # "chr13_1_1000000" : {
            #     "with_a" : ['p', 'a'],
            #     "with_0" : ['p', '0'],
            #     "with_0a" : ['p', '0', 'a']
            # },
            # "chr13_1_20000000" : {
            #     "with_a" : ['p', 'a'],
            #     "with_0" : ['p', '0'],
            #     "with_0a" : ['p', '0', 'a']
            # },
            "chrY_20000000" : {
                "with_Specific" : ['p'],
                "with_a" : ['p', 'a'],
                "with_0" : ['p', '0'],
                "with_0a" : ['p', '0', 'a']
            }
        },
        2:{
            # "chr15_98753195":{
            #     "with_a" : ['m', 'a'],
            #     "with_0" : ['m', '0'],
            #     "with_0a" : ['m', '0', 'a']
            # }
        }
        
    }

    for hap in dict_name_to_triobins.keys():
        for region in dict_name_to_triobins[hap].keys():
            for name, seq_comb in dict_name_to_triobins[hap][region].items():
                dir_out = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/Local_Assembly/ONT_only/Hap{hap}/{region}.{name}"

                path_ont_format = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/Align.Winnowmap2/Region_filtered/Fastq_Filtered/Extract_FastQ.Ref.CHM13.QRY.ONT_UL.Filter_L50kb_QV20_Q30Perc80.triobin_split.triobin_{trio_classify}.{region}.fastq"

                list_path_ont = list(map(lambda trio_classify: path_ont_format.format(region = region, trio_classify = trio_classify), seq_comb))

                assert sum(list(map(lambda path: os.path.exists(path), list_path_ont))) == len(list_path_ont), list_path_ont
                path_ont = ' '.join(list_path_ont)

                curr_time = int(time.time())
                run(
                    path_ont, 
                    dir_out,
                    f"{dir_workspace}/Flye_Local_Assembly/tmp",
                    f"{dir_workspace}/Flye_Local_Assembly/sh",
                    f"{dir_workspace}/Flye_Local_Assembly/log",
                    f"Flye_Local_Assembly.ONT.hap{hap}.{region}.{name}.{curr_time}",
                    32,
                    hostname = [],
                    jobs_wait = []
                )