#%%
import os

import pandas as pd
from Bio.bgzf import BgzfReader


def run(path_vcf, path_depth_on_specific, path_save, min_depth):
    table_vcf, vcf_header = read_vcf_as_table(path_vcf, True)
    table_sv_to_depth = pd.read_csv(path_depth_on_specific, sep = '\t')
    
    dict_chrpos_to_depth = get_chrpos_to_depth(table_sv_to_depth)
    
    table_vcf["depth_on_specific"] = table_vcf.apply(lambda row: dict_chrpos_to_depth[row["#CHROM"]][row["POS"]], axis = 1)
    
    table_vcf_save = table_vcf[table_vcf["depth_on_specific"] < min_depth]
    table_vcf_save = table_vcf_save.drop(columns = ["depth_on_specific"])
    
    save_vcf_for_polishing(table_vcf_save, vcf_header, path_save)

def check_header_len(path_vcf):
    n_header = 0
    header_context = ""
    if path_vcf.endswith(".gz"):
        with BgzfReader(path_vcf, 'rb') as fr:
            for line in fr:
                line = line.decode()
                if line.startswith("#CHROM"):
                    break
                else:
                    header_context += line
                    n_header += 1
    else:
        with open(path_vcf, 'r') as fr:
            for line in fr:
                if line.startswith("#CHROM"):
                    break
                else:
                    header_context += line
                    n_header += 1
    return n_header, header_context

def read_vcf_as_table(path_vcf, return_header = False):
    n_header, header_context = check_header_len(path_vcf)
    table_vcf = pd.read_csv(path_vcf, sep = '\t', skiprows = n_header)

    if return_header:
        return table_vcf, header_context
    else:
        return table_vcf

def get_chrpos_to_depth(table_sv_to_depth_all):
    dict_chrpos_to_depth = dict()
    
    for chrname in table_sv_to_depth_all["#CHROM"].unique():
        dict_chrpos_to_depth[chrname] = dict()
        
        table_sv_to_depth_chrname = table_sv_to_depth_all[table_sv_to_depth_all["#CHROM"] == chrname]
        for pos in table_sv_to_depth_chrname["POS"].unique():
            max_depth = table_sv_to_depth_chrname[table_sv_to_depth_chrname["POS"] == pos]["Depth"].max()
            
            dict_chrpos_to_depth[chrname][pos] = max_depth
    return dict_chrpos_to_depth
            
def save_vcf_for_polishing(table_vcf, header, path_output):
    with open(path_output, 'w') as fw:
        fw.write(header)
        fw.write(get_context_to_write(table_vcf.columns))
        for _, row in table_vcf.iterrows():
            fw.write(get_context_to_write(row.to_list()))

def get_context_to_write(list_values, sep = '\t'):
    return sep.join(list(map(str, list_values))) + '\n'


if __name__ == "__main__":
    dict_hap_to_parental = {
        1:'p',
        2:'m'
    }

    for hap in [1, 2]:
        for tgs in ["ONT", "HiFi"]:
        
            path_vcf = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round8_SNV_Call/SNV_Filtering/Triobinned_for_Evaluation/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.Polish_R8.SV.hap{hap}.{tgs}.triobin_{dict_hap_to_parental[hap]}0a.PASS.NonHomRef.Multiallelic_Norm.VAF05.NonHomRef.vcf.gz"
            
            path_depth_on_specific = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round8_SNV_Call/SNV_Filtering/Triobinned_for_Evaluation/Aligned_Depth.hap{hap}.{tgs}.triobin_{dict_hap_to_parental[hap]}0a.PASS.NonHomRef.Multiallelic_Norm.VAF05.NonHomRef.on_triobin_{dict_hap_to_parental[hap]}.tsv"
            min_depth = 5
            
            path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round8_SNV_Call/SNV_Filtering/Triobinned_for_Evaluation/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.Polish_R8.SV.hap{hap}.{tgs}.triobin_{dict_hap_to_parental[hap]}0a.PASS.NonHomRef.Multiallelic_Norm.VAF05.NonHomRef.Depth{min_depth}onSpecific.vcf"
            
            os.makedirs(os.path.dirname(path_save), exist_ok=True)

            run(path_vcf, path_depth_on_specific, path_save, min_depth)
# %%
if __name__ == "__main__":
    dict_hap_to_parental = {
        1:'p',
        2:'m'
    }
    dict_triobin_exchange = {
        'p':'m',
        'm':'p'
    }

    for hap in [1, 2]:
        path_vcf = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round8_SNV_Call/SNV_Filtering/Triobinned_ONT_HomoAlt/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.Polish_R8.SV.hap{hap}.ONT.triobin_{dict_hap_to_parental[hap]}0a.PASS.NonHomRef.biallelic.VAF08.DP10.Norm.vcf.gz"
        
        path_depth_on_specific = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round8_SNV_Call/SNV_Filtering/Triobinned_ONT_HomoAlt/Aligned_Depth.hap{hap}.ONT.triobin_{dict_hap_to_parental[hap]}0a.PASS.NonHomRef.biallelic.VAF08.DP10.Norm.on_triobin_{dict_hap_to_parental[hap]}.tsv"
        min_depth = 5
        
        path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round8_SNV_Call/SNV_Filtering/Triobinned_ONT_HomoAlt/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.Polish_R8.SV.hap{hap}.ONT.triobin_{dict_hap_to_parental[hap]}0a.PASS.NonHomRef.biallelic.VAF08.DP10.Norm.Depth{min_depth}onSpecific.vcf"

        os.makedirs(os.path.dirname(path_save), exist_ok=True)

        run(path_vcf, path_depth_on_specific, path_save, min_depth)
        