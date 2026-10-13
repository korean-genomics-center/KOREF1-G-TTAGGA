#%%
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).parents[3] / "common"))
from fasta_parser import read_fasta_as_dict, save_dict_fasta_into_file


def run(path_genome, list_additional_genome_dict, path_save):
    dict_fasta = read_fasta_as_dict(path_genome)

    for dict_additional in list_additional_genome_dict:
        dict_fasta_additional = read_fasta_as_dict(dict_additional["Path"])
        for name in dict_additional["Name"]:
            assert dict_fasta.get(name) == None
            dict_fasta[name] = dict_fasta_additional[name]
    
    save_dict_fasta_into_file(dict_fasta, path_save)


if __name__ == "__main__":
    dict_hap_exchange = {
        1:2,
        2:1
    }
    dict_hap_to_sexchr = {
        1:"chrY",
        2:"chrX"
    }

    path_mitochr = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Final_Assembly/Mitochondria_Assembly/final_mitogenome.fasta"
    mitochrname = "chrM"

    path_genome_format = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Final_Assembly/FinalASM.KOREF_S1.Trio_Specific.hap{hap}.Reordered.fasta"

    path_save_format = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/CDS_Fasta/Hap{hap}.Blood_Cellline.CGP_PB/Salmon/Index_Salmon/KOREF1_ASM_Hap{hap}.decoy.fa"

    for hap in [1, 2]:
        anotherhap = dict_hap_exchange[hap]
        run(
            path_genome_format.format(hap = hap), 
            [
                {"Path":path_genome_format.format(hap = anotherhap), "Name":[dict_hap_to_sexchr[anotherhap]]},
                {"Path":path_mitochr, "Name":[mitochrname]}
            ], 
            path_save_format.format(hap = hap))