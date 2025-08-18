#%%
import os
from collections import Counter
from itertools import chain

import pandas as pd
from Bio.bgzf import BgzfReader


#%%
def main(
    path_vcf,
    list_path_vcf_sametrio,
    list_path_vcf_difftrio,
    path_save
):
    table_vcf, vcf_header = read_vcf_as_table(path_vcf, True)
    list_vcfs_sametrio = list(map(read_vcf_as_table, list_path_vcf_sametrio))
    list_vcfs_difftrio = list(map(read_vcf_as_table, list_path_vcf_difftrio))

    table_vcf["VarID"] = get_varids_from_vcf(table_vcf)

    dict_varid_to_gtset_sametrio = get_snvs_with_genotype_combinations_from_vcfs(list_vcfs_sametrio)
    dict_varid_to_gtset_difftrio = get_snvs_with_genotype_combinations_from_vcfs(list_vcfs_difftrio)
    
    # list_varids_observed_both_in_sametrio = list(filter(lambda varid: "./." not in dict_varid_to_gtset_sametrio[varid], dict_varid_to_gtset_sametrio.keys()))
    # list_varids_observed_both_in_difftrio = list(filter(lambda varid: "./." not in dict_varid_to_gtset_difftrio[varid], dict_varid_to_gtset_difftrio.keys()))
    
    table_vcf_polish = evaluate_multiallelic_snvset_with_trio_vcf_snvs(table_vcf.copy(), list(dict_varid_to_gtset_sametrio.keys()), list(dict_varid_to_gtset_difftrio.keys()))
    
    assert set(table_vcf_polish[table_vcf_polish["On_Sametrio"]]["ID"]) == set(table_vcf_polish[table_vcf_polish["On_Difftrio"]]["ID"])
    
    table_vcf_polish_thisasm = table_vcf_polish[table_vcf_polish["On_Sametrio"]]
    table_vcf_polish_thisasm = table_vcf_polish_thisasm.drop(columns = ["VarID", "On_Sametrio", "On_Difftrio", "N_Observed"])
    
    save_vcf_for_polishing(table_vcf_polish_thisasm, vcf_header, path_save)


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

def get_snvs_with_genotype_combinations_from_vcfs(list_vcfs):
    list_dict_varid_to_gt = list()
    total_varid = list()
    for table_vcf in list_vcfs:
        list_varname = get_varids_from_vcf(table_vcf)
        list_genotype = get_genotype_from_vcf(table_vcf, reorder_genotype = True)
        list_dict_varid_to_gt.append(dict(zip(list_varname, list_genotype)))
        total_varid.extend(list_varname)

    dict_varid_to_gt_set = dict()
    for varid in set(total_varid):
        list_genotypes = tuple(map(lambda dict_varid_to_gt: dict_varid_to_gt.get(varid, "./."), list_dict_varid_to_gt))
        dict_varid_to_gt_set[varid] = list_genotypes
    return dict_varid_to_gt_set

def get_varids_from_vcf(table_vcf):
    list_varids = table_vcf.apply(lambda row: f"{row['#CHROM']}:{row['POS']}:{row['REF']}:{row['ALT']}", axis = 1)
    return list_varids

def get_genotype_from_vcf(table_vcf, col_sample = "Sample", reorder_genotype = False):
    list_genotypes = table_vcf[col_sample].apply(lambda val: val.split(':')[0])
    if reorder_genotype:
        list_genotypes = list(map(reorder_genotype_text_with_alphabet_order, list_genotypes))
    return list_genotypes

def reorder_genotype_text_with_alphabet_order(genotype):
    sep = '/' if '/' in genotype else '|'
    list_alleles = list(map(int, genotype.split(sep)))

    return sep.join(list(map(str, sorted(list_alleles))))

def evaluate_multiallelic_snvset_with_trio_vcf_snvs(table_vcf, list_varid_sametrio, list_varid_difftrio):
    set_sametrio_varid = set(list_varid_sametrio)
    set_difftrio_varid = set(list_varid_difftrio)
    
    table_vcf["On_Sametrio"] = table_vcf["VarID"].apply(lambda val: val in set_sametrio_varid)
    table_vcf["On_Difftrio"] = table_vcf["VarID"].apply(lambda val: val in set_difftrio_varid)
    
    table_vcf["N_Observed"] = table_vcf["On_Sametrio"].astype(int)+table_vcf["On_Difftrio"].astype(int)
    table_vcf_singleobserved = table_vcf[table_vcf["N_Observed"] == 1]
    table_vcf_singleobserved["Which_Observed"] = table_vcf_singleobserved.apply(lambda row: "Sametrio" if row["On_Sametrio"] else "Difftrio", axis = 1)
    
    dict_varid_to_observed = dict()
    for _, row in table_vcf_singleobserved.iterrows():
        varid = row["ID"]
        which_obs = row["Which_Observed"]
        
        if dict_varid_to_observed.get(varid) == None:
            dict_varid_to_observed[varid] = list()
        dict_varid_to_observed[varid].append(which_obs)
        
    list_varid_bothobserved = list(filter(lambda varid: len(set(dict_varid_to_observed[varid])) == 2, dict_varid_to_observed.keys()))
    
    table_vcf_polish = table_vcf[table_vcf["ID"].isin(list_varid_bothobserved)]
    return table_vcf_polish

def save_vcf_for_polishing(table_vcf, header, path_output):
    with open(path_output, 'w') as fw:
        fw.write(header)
        fw.write(get_context_to_write(table_vcf.columns))
        for _, row in table_vcf.iterrows():
            fw.write(get_context_to_write(row.to_list()))

def get_context_to_write(list_values, sep = '\t'):
    return sep.join(list(map(str, list_values))) + '\n'
        

#%%
if __name__ == "__main__":
    dict_hap_to_triobin = {
        1:'p',
        2:'m'
    }
    dict_triobin_exchange = {
        'p':'m',
        'm':'p'
    }

    for hap in [1, 2]:
        triobin_same = dict_hap_to_triobin[hap]
        triobin_diff = dict_triobin_exchange[triobin_same]

        path_vcf = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Align_Polished_ASM.From_R1.SV/SNV_Call/SNV_Calling.HiFi.Hap{hap}.Triobin_pm0a/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.R1_Sv.hap{hap}.HiFi.triobin_pm0a.PASS.NonHomRef.multiallelic.GQ25.MinVAF03.SumVAF07.ID_Annot.Norm.vcf.gz"
        path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Align_Polished_ASM.From_R1.SV/SNV_Call/SNV_Calling.HiFi.Hap{hap}.Triobin_pm0a/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.R1_Sv.hap{hap}.HiFi.triobin_pm0a.PASS.NonHomRef.multiallelic.GQ25.MinVAF03.SumVAF07.ID_Annot.Norm.ForPolishing.vcf"

        list_path_vcf_sametrio = [
            f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Align_Polished_ASM.From_R1.SV/SNV_Call/SNV_Calling.HiFi.Hap{hap}.Triobin_{triobin_same}/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.R1_Sv.hap{hap}.HiFi.triobin_{triobin_same}.PASS.NonHomRef.Multiallelic_Norm.VAF05.NonHomRef.vcf.gz",
            f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Align_Polished_ASM.From_R1.SV/SNV_Call/SNV_Calling.ONT.Hap{hap}.Triobin_{triobin_same}/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.R1_Sv.hap{hap}.ONT.triobin_{triobin_same}.PASS.NonHomRef.Multiallelic_Norm.VAF05.NonHomRef.vcf.gz"
        ]
        list_path_vcf_difftrio = [
            f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Align_Polished_ASM.From_R1.SV/SNV_Call/SNV_Calling.HiFi.Hap{hap}.Triobin_{triobin_diff}/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.R1_Sv.hap{hap}.HiFi.triobin_{triobin_diff}.PASS.NonHomRef.Multiallelic_Norm.VAF05.NonHomRef.vcf.gz",
            f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Align_Polished_ASM.From_R1.SV/SNV_Call/SNV_Calling.ONT.Hap{hap}.Triobin_{triobin_diff}/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.R1_Sv.hap{hap}.ONT.triobin_{triobin_diff}.PASS.NonHomRef.Multiallelic_Norm.VAF05.NonHomRef.vcf.gz"
        ]

        main(
            path_vcf,
            list_path_vcf_sametrio,
            list_path_vcf_difftrio,
            path_save
        )
# %%
