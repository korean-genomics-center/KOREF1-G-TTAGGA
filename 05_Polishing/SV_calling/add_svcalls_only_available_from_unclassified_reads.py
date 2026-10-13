#%%
import os

import pandas as pd


def run(path_vcf_without_unclassified, path_vcf_with_unclassified, list_path_depth_of_sv_without_unclassified, path_save, min_depth):
    table_vcf_with_unclassified = read_vcf_as_table(path_vcf_with_unclassified)
    table_sv_to_depth_all = pd.concat(list(map(lambda path: pd.read_csv(path, sep = '\t'), list_path_depth_of_sv_without_unclassified)))
    
    dict_chrpos_to_depth = get_chrpos_to_depth(table_sv_to_depth_all)
    table_vcf_with_unclassified["depth_without_unclassified"] = table_vcf_with_unclassified.apply(lambda row: dict_chrpos_to_depth[row["#CHROM"]][row["POS"]], axis = 1)
    
    table_vcf_with_unclassified_needed = table_vcf_with_unclassified[table_vcf_with_unclassified["depth_without_unclassified"] < min_depth]
    
    list_varids_needed = table_vcf_with_unclassified_needed["ID"].to_list()
    
    merge_vcf_lines_of_specific_varids(path_vcf_without_unclassified, path_vcf_with_unclassified, path_save, list_varids_needed)

def read_vcf_as_table(path_vcf):
    n_header = check_vcf_header_len(path_vcf)
    table_vcf = pd.read_csv(path_vcf, sep = '\t', skiprows = n_header)
    return table_vcf

def check_vcf_header_len(path_vcf):
    n_header = 0

    read_mode = 'rb' if path_vcf.endswith(".gz") else 'r'
    with open(path_vcf, read_mode) as fr:
        for line in fr:
            if hasattr(line, "decode"):
                line = line.decode()
            if line.startswith("#CHROM"):
                break
            else:
                n_header += 1
    return n_header

def get_chrpos_to_depth(table_sv_to_depth_all):
    dict_chrpos_to_depth = dict()
    
    for chrname in table_sv_to_depth_all["#CHROM"].unique():
        dict_chrpos_to_depth[chrname] = dict()
        
        table_sv_to_depth_chrname = table_sv_to_depth_all[table_sv_to_depth_all["#CHROM"] == chrname]
        for pos in table_sv_to_depth_chrname["POS"].unique():
            max_depth = table_sv_to_depth_chrname[table_sv_to_depth_chrname["POS"] == pos]["Depth"].max()
            
            dict_chrpos_to_depth[chrname][pos] = max_depth
    return dict_chrpos_to_depth
            

def merge_vcf_lines_of_specific_varids(path_vcf_original, path_vcf_add, path_save, list_varids_add):
    with open(path_save, 'w') as fw:
        with open(path_vcf_original, 'r') as fr:
            for line in fr:
                fw.write(line)
        
        with open(path_vcf_add, 'r') as fr:
            for line in fr:
                if line.startswith('#'):
                    continue
                else:
                    varid = line.split('\t')[2]
                    if varid in list_varids_add:
                        fw.write(line)


if __name__ == "__main__":
    dict_hap_to_triobin = {
        1:'p',
        2:'m'
    }

    for hap in [1, 2]:

        triobin = dict_hap_to_triobin[hap]
    
        path_vcf_without_unclassified = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round8_SV_Call/SV_Call_without_Unclassified/SV_Calling/Merged_SV/SV_Call.KOREF_PolishedASM_Round8_SV.hap{hap}.hifi_ont.sametrio.filtered.iris.merged.split.homoalt_include.hetalt_only.filt_difftrio_homoref.filt_no_asmsv.vaf_0.5.minsupp_5.call_merged.indel_only.vcf"
        path_vcf_with_unclassified = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round8_SV_Call/SV_Call_with_Unclassified/SV_Calling/Merged_SV/SV_Call.KOREF_PolishedASM_Round8_SV.hap{hap}.hifi_ont.sametrio.filtered.iris.merged.split.homoalt_include.hetalt_only.filt_difftrio_homoref.filt_no_asmsv.vaf_0.5.minsupp_5.call_merged.indel_only.vcf"
        
        list_path_depth_of_sv_without_unclassified = [
            f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round8_SV_Call/SV_Call_Merged/Aligned_Depth_on_TrioSpecific/Aligned_Depth.hap{hap}.ONT_UL.Filter_L10kb_QV15.Trio_Specific_only.triobin_{triobin}.of_SVs_with_Unclassified.tsv",
            f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round8_SV_Call/SV_Call_Merged/Aligned_Depth_on_TrioSpecific/Aligned_Depth.hap{hap}.PacBio_HiFi_ccs_converted.Filter_L1000_QV20.Trio_Specific_only.triobin_{triobin}.of_SVs_with_Unclassified.tsv"
        ]
        min_depth = 5
        
        path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round8_SV_Call/SV_Call_Merged/Merged_Call/SV_Call.KOREF_PolishedASM_Round8_SV.hap{hap}.hifi_ont.sametrio.filtered.iris.final_evaluated.indel_only.additional_support_of_unclassified_reads_on_depth_{min_depth}.vcf"
        
        os.makedirs(os.path.dirname(path_save), exist_ok=True)

        run(path_vcf_without_unclassified, path_vcf_with_unclassified, list_path_depth_of_sv_without_unclassified, path_save, min_depth)
# %%
