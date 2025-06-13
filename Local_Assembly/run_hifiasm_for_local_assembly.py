import pandas as pd
import os, time
import subprocess

PATH_HIFIASM = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/hifiasm_v0.25.0-r726/hifiasm/hifiasm"


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

def run(path_hifi_read, path_save, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list(), **additional_params):    

    cmd = f"""{PATH_HIFIASM} -o {path_save} --n-hap 1 -l0 --telo-m CCCTAA -t{n_jobs}"""
    
    for param_name, param_value in additional_params.items():
        if len(param_name) == 1:
            cmd += f" -{param_name} {param_value}"
        else:
            cmd += f" --{param_name} {param_value}"
    cmd += f" {path_hifi_read}"
        
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
    
    dict_hap_to_region = {
        1:[
            # "chr13:1-1000000",
            # "chr13:1-20000000",
            "chr13:1-5000000"
        ],
        2:[
            # "chr15:94775221"
        ]
    }
    dict_hap_to_parental = {
        1:'p',
        2:'m'
    }
    dict_hap_rev = {
        1:2,
        2:1
    }
    
    for additional_bin in [(),  ('0')]:
        for hap in [1, 2]:
            for region in dict_hap_to_region[hap]:
                additional_name = ''.join(additional_bin)
                if additional_name == '':
                    additional_name = "Specific"

                region_name = region.replace(':', '_').replace('-', '_')
                
                path_hifi_format = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/Align.Winnowmap2/Region_filtered/Fastq_Filtered/Extract_FastQ.Ref.KOREF_DraftASM.hap{dict_hap_rev[hap]}.QRY.PacBio_HiFi_ccs_converted.Filter_L1000_QV20.triobin_split.triobin_{{triobin}}.{region_name}.fastq"
                path_ont_format = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/Align.Winnowmap2/Region_filtered/Fastq_Filtered/Extract_FastQ.Ref.KOREF_DraftASM.hap{dict_hap_rev[hap]}.QRY.ONT_UL.Filter_L50kb_QV20_Q30Perc80.triobin_split.triobin_{{triobin}}.{region_name}.fastq"

                path_hifi = path_hifi_format.format(triobin = dict_hap_to_parental[hap])
                path_ont = path_ont_format.format(triobin = dict_hap_to_parental[hap])
                list_additional_hifi = list(map(lambda triobin: path_hifi_format.format(triobin = triobin), additional_bin))
                list_additional_ont = list(map(lambda triobin: path_ont_format.format(triobin = triobin), additional_bin))

                path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/Local_Assembly/HiFi_ONT/Hap{hap}/{region_name}.Mapped_KOREF_DraftASM_hap{dict_hap_rev[hap]}.with_{additional_name}.l0/KOREF_S1.Local_Assembly.{region_name}.ul_L50kb_QV20_Q30Perc80.triobin_{dict_hap_to_parental[hap]}_{additional_name}.hifi_L1000_QV20.triobin_{dict_hap_to_parental[hap]}_{additional_name}.haploid.l0.asm"

                assert os.path.exists(path_hifi), path_hifi
                assert os.path.exists(path_ont), path_ont
                for addhifi in list_additional_hifi:
                    assert os.path.exists(addhifi), addhifi
                for addont in list_additional_ont:
                    assert os.path.exists(addont), addont

                path_hifi_all = ' '.join([path_hifi]+list_additional_hifi)
                path_ont_hq_all = ','.join([path_ont]+list_additional_ont)

                curr_time = int(time.time())
                run(
                    path_hifi_all,  
                    path_save,
                    f"{dir_workspace}/Hifiasm_Local_Assembly/tmp",
                    f"{dir_workspace}/Hifiasm_Local_Assembly/sh",
                    f"{dir_workspace}/Hifiasm_Local_Assembly/log",
                    f"Hifiasm_Local_Assembly.HiFi_ONT.hap{hap}.l0.with_{additional_name}.{region_name}.Mapped_KOREF_DraftASM.{curr_time}",
                    64,
                    hostname = [],
                    jobs_wait = [],
                    **{
                        "ul" : path_ont_hq_all
                    }
                )