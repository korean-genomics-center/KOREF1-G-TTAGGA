#%%
import subprocess, os, time

PATH_BOWTIEBUILD = "/BiO/Research/Project2/Endangered_Goral/Resources/Tools/Bowtie2/Bowtie2_v2.5.4/bowtie2-2.5.4-linux-x86_64/bowtie2-build"
PATH_SAMTOOLS = "/BiO/Research/Project2/Endangered_Jindo_Dog/Resources/Tools/samtools/samtools_v1.17/samtools-1.17/bin/samtools"
PATH_HICPRO_DIGESTGENOME = "/BiO/Research/Project2/Endangered_Goral/Resources/Tools/HiC-Pro/HiC-Pro_v3.1.0/HiC-Pro/bin/utils/digest_genome.py"


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
    
def copy_file(src, dst):
    subprocess.run(f"cp {src} {dst}", shell = True)
    
def get_config_context(filename_genome, enzymesites_name, list_ligation_motifs, read1_suffix, read2_suffix, n_jobs, ram):
    genome_name = filename_genome.replace(".fasta", '')
    config = f"""# Please change the variable settings below if necessary

#########################################################################
## Paths and Settings  - Do not edit !
#########################################################################

TMP_DIR = tmp
LOGS_DIR = logs
BOWTIE2_OUTPUT_DIR = bowtie_results
MAPC_OUTPUT = hic_results
RAW_DIR = rawdata

#######################################################################
## SYSTEM AND SCHEDULER - Start Editing Here !!
#######################################################################
N_CPU = {n_jobs}
SORT_RAM = {ram}
LOGFILE = hicpro.log

# JOB_NAME = 
# JOB_MEM = 
# JOB_WALLTIME = 
# JOB_QUEUE = 
# JOB_MAIL = 

#########################################################################
## Data
#########################################################################

PAIR1_EXT = {read1_suffix}
PAIR2_EXT = {read2_suffix}

#######################################################################
## Alignment options
#######################################################################

MIN_MAPQ = 10
FORMAT = phred33

BOWTIE2_IDX_PATH = /data/genome/bowtie_index
BOWTIE2_GLOBAL_OPTIONS = --very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder
BOWTIE2_LOCAL_OPTIONS =  --very-sensitive -L 20 --score-min L,-0.6,-0.2 --end-to-end --reorder

#######################################################################
## Annotation files
#######################################################################

REFERENCE_GENOME = {genome_name}
GENOME_SIZE = /data/genome/{filename_genome}.genome

#######################################################################
## Allele specific analysis
#######################################################################

# ALLELE_SPECIFIC_SNP = 

#######################################################################
## Capture Hi-C analysis
#######################################################################

# CAPTURE_TARGET =
# REPORT_CAPTURE_REPORTER = 1

#######################################################################
## Digestion Hi-C
#######################################################################

GENOME_FRAGMENT = /data/genome/{filename_genome}.{enzymesites_name}.bed
LIGATION_SITE = {','.join(list_ligation_motifs)}
# MIN_FRAG_SIZE = 
# MAX_FRAG_SIZE =
# MIN_INSERT_SIZE =
# MAX_INSERT_SIZE =

#######################################################################
## Hi-C processing
#######################################################################

# MIN_CIS_DIST =
GET_ALL_INTERACTION_CLASSES = 1
GET_PROCESS_SAM = 0
RM_SINGLETON = 1
RM_MULTI = 1
RM_DUP = 1

#######################################################################
## Contact Maps
#######################################################################

BIN_SIZE = 20000 40000 150000 500000 1000000
BIN_STEP = 1
MATRIX_FORMAT = upper

#######################################################################
## Normalization
#######################################################################
MAX_ITER = 100
FILTER_LOW_COUNT_PERC = 0.02
FILTER_HIGH_COUNT_PERC = 0
EPS = 0.1
"""
    return config

def run(path_genome, dir_hic_data, list_restriction_motifs, list_ligation_motifs, read1_suffix, read2_suffix, dir_output, dir_tmp, dir_shell, dir_log, jobname, n_jobs, ram, hostname = list(), jobs_wait = list(), **additional_params):
    filename_genome = os.path.basename(path_genome)
    
    enzymesites_name = '_'.join(list(map(lambda val: val.replace('^', ''), list_restriction_motifs)))
    
    # dir_save_results = f"{dir_output}/results"
    dir_save_genome = f"{dir_output}/genome"
    dir_save_genome_index = f"{dir_save_genome}/bowtie_index"
    path_save_config = f"{dir_output}/config-hicpro.txt"
    
    # os.makedirs(dir_save_results, exist_ok=True)
    os.makedirs(dir_save_genome, exist_ok=True)
    os.makedirs(dir_save_genome_index, exist_ok=True)
    
    config_context = get_config_context(filename_genome, enzymesites_name, list_ligation_motifs, read1_suffix, read2_suffix, n_jobs, ram)
    with open(path_save_config, 'w') as fw:
        fw.write(config_context)
    
    cmd = ""
    
    if not os.path.exists(f"{dir_save_genome}/{filename_genome}.genome"):
        cmd += f"""
    cp {path_genome} {dir_save_genome};

    {PATH_SAMTOOLS} faidx {dir_save_genome}/{filename_genome};
    cut -f1,2 {dir_save_genome}/{filename_genome}.fai > {dir_save_genome}/{filename_genome}.genome;
    """
    
    filename_genome_with_suffix = filename_genome.replace(".fasta", '')
    if not os.path.exists(f"{dir_save_genome}/bowtie_index/{filename_genome_with_suffix}.rev.2.bt2"):
        cmd += f"""
{PATH_BOWTIEBUILD} --threads {n_jobs} {dir_save_genome}/{filename_genome} {dir_save_genome}/bowtie_index/{filename_genome_with_suffix};

"""

    if not os.path.exists(f"{dir_save_genome}/{filename_genome}.{enzymesites_name}.bed"):
        cmd += f"""
python {PATH_HICPRO_DIGESTGENOME} -r {' '.join(list_restriction_motifs)} -o {dir_save_genome}/{filename_genome}.{enzymesites_name}.bed {dir_save_genome}/{filename_genome}
"""

    cmd_hicpro = f"""
docker run --ipc=host \\
    -v {dir_hic_data}:/input/sample \\
    -v {dir_output}:/data \\
    nservant/hicpro:latest \\
    /HiC-Pro_3.0.0/bin/HiC-Pro \\
    -i /input \\
    -o /data/results \\
    -c /data/config-hicpro.txt \\
"""

    cmd += cmd_hicpro
    
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

if __name__ == "__main2__":
    dir_workspace = "/BiO/Research/Project2/Endangered_Goral/Analysis"
    
    list_restriction_motifs = ["^GATC", "G^ANTC", "C^TNAG", "T^TAA"]
    list_ligation_motifs = ["GATCGATC","GATCANTC","GATCTNAG","GATCTAA","GANTGATC","GANTANTC","GANTTNAG","GANTTAA","CTNAGATC","CTNAANTC","CTNATNAG","CTNATAA","TTAGATC","TTAANTC","TTATNAG","TTATAA"]
    
    for hap in [1]:
        path_genome = f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Initial_Assembly/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.asm.hic.hap{hap}.p_ctg.fasta"
        
        dir_hic_data = "/BiO/Research/Project2/Endangered_Goral/Results/Input_Data_Filtering/ShortRead/HiC/KoGoral/Merged_Fastq"
        
        read1_suffix = "R1"
        read2_suffix = "R2"
        
        dir_output = f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/YaHS_Scaffolding/HiC_Alignment.4Enz/Hap{hap}"
        
        curr_time = int(time.time())
        run(
            path_genome,
            dir_hic_data,
            list_restriction_motifs,
            list_ligation_motifs,
            read1_suffix,
            read2_suffix,
            dir_output,
            f"{dir_workspace}/HiC_Alignment/tmp",
            f"{dir_workspace}/HiC_Alignment/sh",
            f"{dir_workspace}/HiC_Alignment/log",
            f"HiC_Alignment.4Enz.HiCPro.KoGoral.DraftASM.Hap{hap}.{curr_time}",
            n_jobs = 192,
            ram = "4000M",
            hostname = [],
            jobs_wait = []
        )
#%%

if __name__ == "__main__":
    dir_workspace = "/BiO/Research/Project2/Endangered_Goral/Analysis"
    
    list_restriction_motifs = ["^GATC", "G^ANTC", "C^TNAG", "T^TAA"]
    list_ligation_motifs = ["GATCGATC","GATCANTC","GATCTNAG","GATCTAA","GANTGATC","GANTANTC","GANTTNAG","GANTTAA","CTNAGATC","CTNAANTC","CTNATNAG","CTNATAA","TTAGATC","TTAANTC","TTATNAG","TTATAA"]
    
    for hap in [1]:
        path_genome = f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/HiC_Scaffolding_with_additional_contigs/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.asm.hic.hap{hap}.p_ctg.with_additional_teloseq.with_unclustered_rawONT.fasta"
        
        dir_hic_data = "/BiO/Research/Project2/Endangered_Goral/Results/Input_Data_Filtering/ShortRead/HiC/KoGoral/Merged_Fastq"
        
        read1_suffix = "R1"
        read2_suffix = "R2"
        
        dir_output = f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/HiC_Scaffolding_with_additional_contigs/HiC_Alignment.4Enz.with_raw_unclustered_read"
        
        curr_time = int(time.time())
        run(
            path_genome,
            dir_hic_data,
            list_restriction_motifs,
            list_ligation_motifs,
            read1_suffix,
            read2_suffix,
            dir_output,
            f"{dir_workspace}/HiC_Alignment/tmp",
            f"{dir_workspace}/HiC_Alignment/sh",
            f"{dir_workspace}/HiC_Alignment/log",
            f"HiC_Alignment.4Enz.HiCPro.KoGoral.DraftASM.AdditionalTelo.with_raw_unclustered_read.Hap{hap}.{curr_time}",
            n_jobs = 256,
            ram = "4000M",
            hostname = [],
            jobs_wait = []
        )