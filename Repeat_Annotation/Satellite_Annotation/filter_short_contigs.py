#%%
import sys

sys.path.append("/BiO/Research/Project2/Endangered_Jindo_Dog/Resources/Scripts/KOGIC_KOREF_AGGA")
from fasta_parser import read_fasta_as_dict, save_specific_seqnames_from_fasta_dict


# %%
def run(path_srf_contig, min_length, path_save):
    dict_srf_contig = read_fasta_as_dict(path_srf_contig)
    
    list_contigs_over_min_length = list(filter(lambda name: len(dict_srf_contig[name]) >= min_length, dict_srf_contig.keys()))
    
    print(list_contigs_over_min_length)
    
    save_specific_seqnames_from_fasta_dict(dict_srf_contig, list_contigs_over_min_length, path_save)
    

if __name__ == "__main__":
    for hap in [1, 2]:
        for ci in [10, 20]:
            path_srf_contig = f"/BiO/Research/Project2/Endangered_Jindo_Dog/Results/Repeat_Annotation/Satellite_Repeat/run_srf/Hap{hap}/Satellite_repeat_detection.J495799_Child.Hap{hap}.k101.ci{ci}.srf_output.fa"
            
            min_length = 100
            path_save = f"/BiO/Research/Project2/Endangered_Jindo_Dog/Results/Repeat_Annotation/Satellite_Repeat/run_srf/Hap{hap}/Satellite_repeat_detection.J495799_Child.Hap{hap}.k101.ci{ci}.srf_output.min_length_{min_length}.fa"
            
            run(path_srf_contig, min_length, path_save)
        
if __name__ == "__main2__":
    for k in [101, 151]:
        for ci in [10, 20]:
            path_srf_contig = f"/BiO/Research/Project2/Endangered_Jindo_Dog/Results/Repeat_Annotation/Satellite_Repeat/run_srf/Mixed/Satellite_repeat_detection.J495799_Child.Hap_Mixed.k{k}.ci{ci}.srf_output.fa"
            
            min_length = 100
            path_save = f"/BiO/Research/Project2/Endangered_Jindo_Dog/Results/Repeat_Annotation/Satellite_Repeat/run_srf/Mixed/Satellite_repeat_detection.J495799_Child.Hap_Mixed.k{k}.ci{ci}.srf_output.min_length_{min_length}.fa"
            
            run(path_srf_contig, min_length, path_save)
# %%
