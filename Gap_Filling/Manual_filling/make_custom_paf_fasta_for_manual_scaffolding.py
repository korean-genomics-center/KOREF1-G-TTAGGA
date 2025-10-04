#%%
import sys,os
from pathlib import Path

sys.path.append(str(Path(__file__).parents[2]))
from fasta_parser import read_fasta_as_dict, save_dict_fasta_into_file


def run(
    dict_name_to_contignames,
    dict_name_to_contigpath,
    path_paf_format,
    path_save_prefix
):
    os.makedirs(os.path.dirname(path_save_prefix), exist_ok=True)
    
    path_save_fasta = f"{path_save_prefix}.contigs_gathered.fasta"
    path_save_paf = f"{path_save_prefix}.contigs_gathered.aligned.paf"
    
    dict_new_fasta = dict()
    for name, path_fasta in dict_name_to_contigpath.items():
        dict_fasta_original = read_fasta_as_dict(path_fasta)
        for contigname in dict_name_to_contignames[name]:
            new_contigname = f"{name}__{contigname}"
            dict_new_fasta[new_contigname] = dict_fasta_original[contigname]
    save_dict_fasta_into_file(dict_new_fasta, path_save_fasta)
    
    with open(path_save_paf, 'w') as fw:
        for name, list_contignames in dict_name_to_contignames.items():
            path_paf_original = path_paf_format.format(name = name)
            with open(path_paf_original, 'r') as fr:
                for line in fr:
                    contigname_on_paf = line.split('\t')[0]
                    if contigname_on_paf in list_contignames:
                        fw.write(line.replace(contigname_on_paf, f"{name}__{contigname_on_paf}"))
                


#%%
if __name__ == "__main__":
    path_paf_format = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon_wo_RawONT/Round8.Manual_Scaff/Hap1/Prealignment/Align.Original_Contig/Align.REF_KOREFDraftASM_Hap1_Split.QRY_{name}.maphifiPw5s100c.filt_match100_idy0.2.paf"
    
    dict_name_to_contigpath = {
        "Flye_Binned_ONT":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/ONT_UL.Filter_L50kb_QV20_Q30Perc80/Hap1/assembly.fasta",
        "Flye_Binned_HiFi_with0" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20.with_0/Hap1/assembly.fasta"
    }
    
    dict_name_to_contignames = {
        "Flye_Binned_ONT" : [
            "contig_603",
            "contig_1389"
        ],
        "Flye_Binned_HiFi_with0":[
            "contig_11070"
        ]
    }
    
    path_save_prefix = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon_wo_RawONT/Round8.Manual_Scaff/Hap1/TGSGapcloser_with_custom_paf/custom_data_for_gapclosing.maphifiPw5s100c"
    
    run(
        dict_name_to_contignames,
        dict_name_to_contigpath,
        path_paf_format,
        path_save_prefix
    )