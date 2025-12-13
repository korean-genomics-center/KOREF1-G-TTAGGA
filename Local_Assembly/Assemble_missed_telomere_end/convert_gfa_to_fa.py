#%%
import subprocess


list_contignames = [
    "h1tg000011l",
    "h1tg000012l",
    "h1tg000026l"
]

path_gfa_foramt = "/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Telomere_from_RawRead/Assembly_Elongation/from_round1/{contig}_elongate_1_end.Hifiasm.HiFi/{contig}_elongate_1_end.hifi.haploid.l0.asm.bp.p_ctg.gfa"

for contig in list_contignames:
    path_gfa = path_gfa_foramt.format(contig = contig)
    path_save = path_gfa.replace(".gfa", ".fasta")
    
    cmd = f"""awk '/^S/{{print ">"$2"\\n"$3}}' {path_gfa} > {path_save} """
    
    subprocess.run(cmd, shell = True)
    print(cmd)