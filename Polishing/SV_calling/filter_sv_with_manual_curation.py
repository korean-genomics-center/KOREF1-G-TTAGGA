#%%
import json
import os


def run(path_vcf, path_curation, path_save):
    dict_curation = read_json(path_curation)
    
    list_varids_filtout = list(filter(lambda name: dict_curation[name] == False, dict_curation.keys()))
    save_vcf_without_specific_varids(path_vcf, path_save, list_varids_filtout)

def read_json(path_json):
    with open(path_json, 'r') as fr:
        dict_val = json.load(fr)
    return dict_val

def save_vcf_without_specific_varids(path_in, path_out, list_varids_filtout):
    list_filted = list()
    with open(path_in, 'r') as fr, open(path_out, 'w') as fw:
        for line in fr:
            if line.startswith('#'):
                fw.write(line)
            else:
                varid = line.split('\t')[2]
                if varid in list_varids_filtout:
                    list_filted.append(varid)
                    continue
                else:
                    fw.write(line)
    
    assert len(set(list_filted)) == len(list_varids_filtout)
    

if __name__ == "__main__":
    hap = 2
    
    path_vcf = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round4_SV_Call/SV_Call_Merged/Merged_Call/SV_Call.KOREF_PolishedASM_Round4_SV.hap{hap}.hifi_ont.sametrio.filtered.iris.final_evaluated.indel_only.additional_support_of_unclassified_reads_on_depth_5.vcf"
    path_curation = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round4_SV_Call/SV_Call_Merged/Merged_Call/manual_curation_large_indel.hap{hap}.json"
    path_save =  f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round4_SV_Call/SV_Call_Merged/Merged_Call/SV_Call.KOREF_PolishedASM_Round4_SV.hap{hap}.hifi_ont.sametrio.filtered.iris.final_evaluated.indel_only.additional_support_of_unclassified_reads_on_depth_5.manual_curation.vcf"
    
    run(path_vcf, path_curation, path_save)