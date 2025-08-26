#%%
import subprocess, os, time

PATH_DORADO = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/dorado/dorado-0.9.6-linux-x64/bin/dorado"

def get_context_for_qsub_submission(cmd, jobname, stderr_path, stdout_path, nthreads, hold_jids = list()):
    #Header
    script_header = f"#!bin/sh\n"
    script_header += f"#$ -N {jobname}\n"
    script_header += f"#$ -S /bin/bash\n"
    script_header += f"#$ -V\n"
    if len(hold_jids) > 0:
        script_header += f"#$ -hold_jid {','.join(hold_jids)}\n"
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
    
def run_job(path_fastq, path_save, sequencing_kit, jobsuffix, dir_tmp, dir_shell, dir_log, n_jobs = 50):
    cmd = """{PATH_DORADO} trim --sequencing-kit {sequencing_kit} --threads {n_jobs} --emit-fastq {path_fastq} > {path_save}""".format(
        PATH_DORADO = PATH_DORADO,
        path_fastq = path_fastq,
        path_save = path_save,
        sequencing_kit = sequencing_kit,
        n_jobs = n_jobs
    )
    
    jobname = f"Trim_Adapter_ONT_{jobsuffix}.{int(time.time())}"
    
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
            list()
        )
    )


if __name__ == "__main2__":
    workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis/Trim_ONT_Adapter"
    
    path_fastq = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/RawData_WGS/KOREF-S1/Blood/ONT/20240517_1158_P2S-01506-B_PAW01562_215af4e8.fastq"
    path_save = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/ONT/Adapter_Filtered/20240517_1158_P2S-01506-B_PAW01562_215af4e8.Adapt_Filt.fastq"
    
    sequencing_kit = "SQK-ULK114"
    
    run_job(
        path_fastq,
        path_save,
        sequencing_kit,
        "KOREF_ONT_20240517_1158_P2S-01506-B_PAW01562_215af4e8",
        os.path.join(workspace, "tmp"),
        os.path.join(workspace, "sh"),
        os.path.join(workspace, "log"),
        8
    )


if __name__ == "__main__":
    list_seq_name = [
        "20240629_1736_P2S-01506-A_PAW01687_bf35a486",
        "20240629_1736_P2S-01506-B_PAU99972_b4ed6da3",
        "20240811_1107_P2S-01506-A_PAW77551_0cd923fa",
        "20240816_0603_P2S-01501-B_PAW77472_d96fb6c9",
        "20240817_0645_P2S-01501-B_PAW77472_c8fb51b0",
        "20240818_1702_P2S-01506-B_PAW88380_5c7050f8",
        "20240920_1805_P2S-01506-A_PAY14781_7526b2f4",
        "20240923_1616_P2S-01506-A_PAY08433_74256187",
        "20240923_1641_P2S-01506-B_PAY03655_fb6f2a12",
        "20240924_2038_P2S-01506-A_PAY08433_7a23907e",
        "20240924_2038_P2S-01506-B_PAY03655_f02a3d27",
        "20240926_1205_P2S-01506-A_PAY09654_90f482e1",
        "20240926_1205_P2S-01506-B_PAY09817_4744eb86",
        "20240930_1539_P2S-01506-B_PAY10807_076f0b43",
        "20241002_1732_P2S-01506-A_PAY10807_3ab61abe",
        "20241011_1754_P2S-01506-A_PAW99288_c001806f",
        "20241012_1341_P2S-01506-A_PAW99288_c255cb62",
        "20241012_1341_P2S-01506-B_PAY01199_fef5a120",
        "20241021_1517_P2S-01506-A_PAW99288_e821305c",
        "20241021_1517_P2S-01506-B_PAY01199_b712d259",
        "20241022_1541_P2S-01506-A_PAW99013_6c642a53",
        "20241022_1541_P2S-01506-B_PAY07547_1caa91a9",
        "20241031_1212_P2S-01506-A_PAY83904_baf56bf4",
        "20241031_1220_P2S-01506-B_PAY74072_775e6f78",
        "20241031_1440_P2S-01503-A_PAY74208_5cb0c269",
        "20241031_1440_P2S-01503-B_PAY74099_bed3405c",
        "20241101_1326_P2S-01506-A_PAY83904_5bb2fcbf",
        "20241101_1326_P2S-01506-B_PAY74072_95afdfdf",
        # "20240510_1615_P2S-01506-A_PAU90966_ad72ac81_duplex",
        # "20240510_1615_P2S-01506-A_PAU90966_ad72ac81_simplex",
        # "20240510_1615_P2S-01506-B_PAU91096_39bb4496_duplex",
        # "20240510_1615_P2S-01506-B_PAU91096_39bb4496_simplex",
        # "20240603_1341_P2S-01506-A_PAU69619_32394ac4_duplex",
        # "20240603_1341_P2S-01506-A_PAU69619_32394ac4_simplex",
        # "20240603_1341_P2S-01506-B_PAU69571_419b8c25_duplex",
        # "20240603_1341_P2S-01506-B_PAU69571_419b8c25_simplex",
        # "20240603_1553_P2S-01506-A_PAU69619_a90db93e_duplex",
        # "20240603_1553_P2S-01506-A_PAU69619_a90db93e_simplex",
        # "20240603_1553_P2S-01506-B_PAU69571_1464b73b_duplex",
        # "20240603_1553_P2S-01506-B_PAU69571_1464b73b_simplex",
        # "20240603_1925_P2S-01506-A_PAU69619_b2756a99_duplex",
        # "20240603_1925_P2S-01506-A_PAU69619_b2756a99_simplex",
        # "20240603_1925_P2S-01506-B_PAU69571_bd44489b_duplex",
        # "20240603_1925_P2S-01506-B_PAU69571_bd44489b_simplex",
        # "20240603_2030_P2S-01506-A_PAU69619_fad37a7f_duplex",
        # "20240603_2030_P2S-01506-A_PAU69619_fad37a7f_simplex",
        # "20240603_2030_P2S-01506-B_PAU69571_54256320_duplex",
        # "20240603_2030_P2S-01506-B_PAU69571_54256320_simplex",
        # "20240604_1338_P2S-01506-A_PAU90972_437fd6b6_duplex",
        # "20240604_1338_P2S-01506-A_PAU90972_437fd6b6_simplex",
        # "20240604_1338_P2S-01506-B_PAU91173_a1c4c20d_duplex",
        # "20240604_1338_P2S-01506-B_PAU91173_a1c4c20d_simplex",
        # "20240606_1426_P2S-01506-A_PAU91084_50a9eb35_duplex",
        # "20240606_1426_P2S-01506-A_PAU91084_50a9eb35_simplex",
        # "20240606_1426_P2S-01506-B_PAU91079_79717851_duplex",
        # "20240606_1426_P2S-01506-B_PAU91079_79717851_simplex",
        # "20240608_1349_P2S-01506-A_PAU91135_57288ce9_duplex",
        # "20240608_1349_P2S-01506-A_PAU91135_57288ce9_simplex",
        # "20240608_1349_P2S-01506-B_PAU68508_19aa0ef8_duplex",
        # "20240608_1349_P2S-01506-B_PAU68508_19aa0ef8_simplex",
        # "20240611_1318_P2S-01506-A_PAU91094_93be457b_duplex",
        # "20240611_1318_P2S-01506-A_PAU91094_93be457b_simplex"
    ]
    
    workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis/Trim_ONT_Adapter"
    # path_fastq = "/BiO/Access/yoonsung/Research/Test_m6A_Methylation/Resources/RawData/KOREF-S1_Blood_WTS_Sequel-II_20230519_FC00.subreads.bam"
    # path_save = "/BiO/Access/yoonsung/Research/Test_m6A_Methylation/Results/Subread_to_CCS/KOREF-S1_Blood_WTS_Sequel-II_20230519_FC00.ccs.bam"
    dir_save = "/BiO/Access/yoonsung/Research/KOREF_Assembly/Results/PacBio/Subread_to_CCS"
    
    for seqname in list_seq_name:
        path_fastq = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/RawData_WGS/KOREF-S1/Blood/ONT/{seqname}.fastq"
        path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/ONT/Adapter_Filtered/{seqname}.Adapt_Filt.fastq"
        
        sequencing_kit = "SQK-ULK114"
        
        run_job(
            path_fastq,
            path_save,
            sequencing_kit,
            seqname,
            os.path.join(workspace, "tmp"),
            os.path.join(workspace, "sh"),
            os.path.join(workspace, "log"),
            8
        )
#%%