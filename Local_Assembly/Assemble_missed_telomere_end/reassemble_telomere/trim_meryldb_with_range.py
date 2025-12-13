#%%
import subprocess

PATH_MERYL = "/BiO/Research/Project2/Endangered_Jindo_Dog/Resources/Tools/meryl/meryl_v1.4.1/meryl-1.4.1/bin/meryl"

ratio = 65
k = 17

min_cnt = 12
max_cnt = 42

path_db = f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/Classify_ONT_reads/Meryl.unaligned_Telorepeat_ONT_reads_ratio{ratio}.k{k}"

path_out1 = f"{path_db}.min_{min_cnt}"
path_out2 = f"{path_db}.min_{min_cnt}.max_{max_cnt}"

cmd = f"""
{PATH_MERYL} greater-than {min_cnt-1} output {path_out1} {path_db};
{PATH_MERYL} less-than {max_cnt+1} output {path_out2} {path_out1};
"""

subprocess.run(cmd, shell = True)


#%%
cmd_print_kmer = f"{PATH_MERYL} print {path_out2} > {path_out2}.print"
subprocess.run(cmd_print_kmer, shell = True)