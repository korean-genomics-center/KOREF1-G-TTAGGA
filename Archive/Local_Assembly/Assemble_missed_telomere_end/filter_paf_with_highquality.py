#%%
from joblib import Parallel, delayed
import subprocess

minidy = 0.90
minalign = 1000
maxhang = 1000
mina2h = 0.90

path_script = "/BiO/Research/Project2/Endangered_Jindo_Dog/Resources/Tools/TGS-Gapcloser/TGS-Gapcloser_fork/TGS-GapCloser/filter_hanged_candidate.py"


path_paf_format = "/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/From_ONT/RawReads_Aligned_ONT/from_round1/Minimap2_Align/Align_Minimap2.Ref_{contig}.R1.QRY_Raw_ONT.mapontN20cw5U101000000e0m100r2k.endaligned1000.tmp.paf"
list_contigs = [
    "h1tg000012l",
    # "h1tg000026l"
]
list_path_paf = list(map(lambda name: path_paf_format.format(contig = name), list_contigs))


def run_filtering_script(path_paf):
    path_out = f"{path_paf.replace('.paf', '')}.idy{int(minidy*100)}align{minalign}hang{maxhang}a2h{int(mina2h*100)}.paf"
    
    cmd = f"python {path_script} --paf {path_paf} --out {path_out} --idy {minidy} --align {minalign} --hang_wrong {maxhang} --align_ratio {mina2h}"
    subprocess.run(cmd, shell = True)
    print(f"Finished: {path_paf}", flush = True)


with Parallel(n_jobs = len(list_contigs)) as parallel:
    parallel(delayed(run_filtering_script)(path_paf_contig) for path_paf_contig in list_path_paf)