#%%
from joblib import Parallel, delayed
import subprocess,os

minidy = 0.95
minalign = 1000
maxhang = 10000
mina2h = 0.95

minreadlength = 100000
minalignlength = 75000

path_script = "/BiO/Research/Project2/Endangered_Jindo_Dog/Resources/Tools/TGS-Gapcloser/TGS-Gapcloser_fork/TGS-GapCloser/filter_hanged_candidate.py"

#%%
def run_filtering_script(path_paf):
    path_out = f"{path_paf.replace('.paf', '')}.idy{int(minidy*100)}align{minalign}hang{maxhang}a2h{int(mina2h*100)}.paf"
    
    cmd = f"python {path_script} --paf {path_paf} --out {path_out} --idy {minidy} --align {minalign} --hang_wrong {maxhang} --align_ratio {mina2h}"
    subprocess.run(cmd, shell = True)
    print(f"Finished: {path_paf}", flush = True)
    
#%%
for alignopt in ["a2h80idy95"]:
    for ratio in [60]:
        
        print("Start: ", alignopt, ratio)
        
        # path_paf = f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Telomere_with_read_classification/Cluster_Unaligned_ONT_Reads_by_Alignment/Minimap2_Align/Align_Minimap2.Telomeric_Read_Self_bothhap_{alignopt}_ratio_{ratio}.mapontPcw5U11000000e0.paf"
        path_paf = f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Initial_Assembly/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_Detection/Fasta_Split_by_Telomere/Minimap2_Align/Align_Minimap2.Ref_nonT2TcontigHap2.Qry_nonT2TcontigHap2.maphifiPcw5s100.paf"



        run_filtering_script(path_paf)

        cmd_filter_long_align = f"""
        awk '$2>{minreadlength} && $7>{minreadlength} {{print}}' {path_paf.replace('.paf', '')}.idy{int(minidy*100)}align{minalign}hang{maxhang}a2h{int(mina2h*100)}.paf > {path_paf.replace('.paf', '')}.idy{int(minidy*100)}align{minalign}hang{maxhang}a2h{int(mina2h*100)}.lenover_{minreadlength}.paf;
        awk '$11>{minalignlength} {{print}}' {path_paf.replace('.paf', '')}.idy{int(minidy*100)}align{minalign}hang{maxhang}a2h{int(mina2h*100)}.lenover_{minreadlength}.paf > {path_paf.replace('.paf', '')}.idy{int(minidy*100)}align{minalign}hang{maxhang}a2h{int(mina2h*100)}.lenover_{minreadlength}.alignover_{minalignlength}.paf
        """

        # subprocess.run(cmd_filter_long_align, shell = True)
        
        print(os.path.basename(path_paf))