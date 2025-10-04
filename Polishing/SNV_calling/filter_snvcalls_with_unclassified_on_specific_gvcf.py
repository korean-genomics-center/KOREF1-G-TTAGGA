#%%
import os, math
from collections import Counter
from itertools import chain

from joblib import Parallel, delayed
import pandas as pd
from Bio.bgzf import BgzfReader


#%%
def main(
    path_vcf,
    list_path_gvcfs,
    path_save,
    path_save_stat,
    n_jobs = 30
):
    table_vcf, vcf_header = read_vcf_as_table(path_vcf, True)
    
    table_vcf["VarID"] = get_varids_from_vcf(table_vcf)
    
    list_varids_na_on_gvcf = filter_na_varids_on_specific_only_gvcf(table_vcf["VarID"].to_list(), list_path_gvcfs, n_jobs)
    
    table_vcf_save = table_vcf[table_vcf["VarID"].isin(list_varids_na_on_gvcf)]
    table_vcf_save = table_vcf_save.drop(columns = ["VarID"])
    
    with open(path_save_stat, 'w') as fw:
        fw.write(f"Total SNVs: {table_vcf.shape[0]}\n")
        fw.write(f"SNVs NA on Specific-only gVCF: {len(list_varids_na_on_gvcf)}\n")
        
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
    
def get_varids_from_vcf(table_vcf):
    list_varids = table_vcf.apply(lambda row: f"{row['#CHROM']}:{row['POS']}:{row['REF']}:{row['ALT']}", axis = 1)
    return list_varids

def filter_na_varids_on_specific_only_gvcf(list_varids_check, list_path_gvcf, n_jobs):
    dict_varids_to_gts_on_gvcf = check_varids_on_gvcf(list_varids_check, list_path_gvcf, n_jobs)

    list_varids_na_on_gvcf = list()
    for varid, list_gts_on_gvcf in dict_varids_to_gts_on_gvcf.items():
        set_gts_on_gvcf = set(list_gts_on_gvcf)
        if set_gts_on_gvcf <= {"./.", ".|."}:
            list_varids_na_on_gvcf.append(varid)
    return list_varids_na_on_gvcf

def check_varids_on_gvcf(list_varids_check, list_path_gvcf, n_jobs = 2):

    list_varids_batched_by_chrname = split_varids_by_chrnames_with_sorting(list_varids_check)

    with Parallel(n_jobs = min(len(list_varids_batched_by_chrname), n_jobs)) as parallel:
        list_dict_varid_to_gvcfcall = parallel(delayed(check_varids_on_gvcf_single_chromosome)(
            list_varids_check_batch, list_path_gvcf
        )for list_varids_check_batch in list_varids_batched_by_chrname)

    dict_varid_to_gvcfcall = dict()
    for dict_varid_to_gvcfcall_batch in list_dict_varid_to_gvcfcall:
        assert len(set(dict_varid_to_gvcfcall.keys()) & set(dict_varid_to_gvcfcall_batch.keys())) == 0
        dict_varid_to_gvcfcall.update(dict_varid_to_gvcfcall_batch)
    
    return dict_varid_to_gvcfcall

def split_varids_by_chrnames_with_sorting(list_varids):
    list_chrnames_on_varids = list(set(list(map(lambda val: val.split(':')[0], list_varids))))

    list_varid_batches_by_chrname = list()
    for chrname in list_chrnames_on_varids:
        list_varids_of_chrname = list(filter(lambda varid: varid.split(':')[0] == chrname, list_varids))
        list_varids_of_chrname_sorted = sorted(list_varids_of_chrname, key = lambda val: int(val.split(':')[1]))
        list_varid_batches_by_chrname.append(list_varids_of_chrname_sorted)
    return list_varid_batches_by_chrname

def check_varids_on_gvcf_single_chromosome(list_varids_check_batch, list_path_gvcf):
    chrname_of_this_batch = list_varids_check_batch[0].split(':')[0]

    list_gvcfs_region_difftrio = list(map(read_gvcf_region, list_path_gvcf))
    for dict_gvcf_region in list_gvcfs_region_difftrio:
        verify_gvcf_regions_sorted_of_regions_on_single_chromosome(chrname_of_this_batch, list(dict_gvcf_region[chrname_of_this_batch]))
    
    dict_varid_to_gvcfcall = dict()

    dict_gvcf_ind_to_start_searching = dict()
    for ind in range(len(list_gvcfs_region_difftrio)):
        dict_gvcf_ind_to_start_searching[ind] = 0

    for varid in list_varids_check_batch:
        dict_varid_to_gvcfcall[varid] = list()
        varid_chrname, varid_pos, _, _ = varid.split(':')
        varid_pos = int(varid_pos)

        for ind_gvcf, dict_gvcf_region in enumerate(list_gvcfs_region_difftrio):
            list_regions_on_gvcf = list(dict_gvcf_region[varid_chrname].keys())
            for ind_region, region in enumerate(list_regions_on_gvcf[dict_gvcf_ind_to_start_searching[ind_gvcf]:], start = dict_gvcf_ind_to_start_searching[ind_gvcf]):
                if len(region) == 0:
                    raise Exception(f"{varid_chrname}")
                if check_is_value_in_range(region, varid_pos):
                    dict_gvcf_ind_to_start_searching[ind_gvcf] = ind_region
                    break
                elif region[0] > varid_pos:
                    raise Exception(f"{varid} not found in gvcf regions: Index {ind_gvcf}")

            ind_region_containing_position = ind_region

            list_ranges_ind_check = range(max(0, ind_region_containing_position-100), min(ind_region_containing_position+101, len(list_regions_on_gvcf)-1)+1)
            list_regions_check = list(map(lambda ind: list_regions_on_gvcf[ind], list_ranges_ind_check))
            list_regions_containing_pos = list(filter(lambda region: check_is_value_in_range(region, varid_pos), list_regions_check))

            if ind_region_containing_position-100 > 0:
                assert list_regions_check[0] != list_regions_containing_pos[0], list_regions_check
            if ind_region_containing_position+101 < len(list_regions_on_gvcf)-1:
                assert list_regions_check[-1] != list_regions_containing_pos[-1], list_regions_check

            list_gts_containing_pos = list(map(dict_gvcf_region[varid_chrname].__getitem__, list_regions_containing_pos))

            assert len(list_gts_containing_pos) > 0

            dict_varid_to_gvcfcall[varid].extend(list_gts_containing_pos)
    print(chrname_of_this_batch, "Finished", flush = True)
    return dict_varid_to_gvcfcall


def read_gvcf_region(path_vcf):
    if path_vcf.endswith(".gz"):
        filereader = BgzfReader
        mode = "rb"
        decode = True
    else:
        filereader = open
        mode = "r"
        decode = False

    dict_chrname_to_gvcf_calls = dict()

    with filereader(path_vcf, mode) as fr:
        for line in fr:
            if decode: line = line.decode()
            if line.startswith("#CHROM"):
                break
        for line in fr:
            if decode: line = line.decode()
            list_values = line.split('\t')
            chrname = list_values[0]
            ref = list_values[3]
            start = int(list_values[1])
            end = list_values[7].split('=')[-1]
            end = start+len(ref)-1 if end == '.' else int(end)
            gt = list_values[-1].split(':')[0]

            if dict_chrname_to_gvcf_calls.get(chrname) == None:
                dict_chrname_to_gvcf_calls[chrname] = dict()

            poskey = (start, end)
            assert dict_chrname_to_gvcf_calls[chrname].get(poskey) == None, f"{chrname}: {poskey}"
            dict_chrname_to_gvcf_calls[chrname][poskey] = gt
    return dict_chrname_to_gvcf_calls

def verify_gvcf_regions_sorted_of_regions_on_single_chromosome(chrname, list_regions):
    region_before = list_regions[0]
    for region in list_regions[1:]:
        assert region[0] > region_before[0], f"{chrname}, {region_before}, {region}"
        region_before = region
        
def check_is_value_in_range(ranges, value):
    return value >= ranges[0] and value <= ranges[1]

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
    dict_hap_to_parental = {
        1:'p',
        2:'m'
    }
    dict_triobin_exchange = {
        'p':'m',
        'm':'p'
    }


    for hap in [1, 2]:
        path_vcf = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round3_SNV_Call/SNV_Filtering/Triobinned_ONT_HomoAlt/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.Polish_R3.SV.hap{hap}.ONT.triobin_{dict_hap_to_parental[hap]}0a.PASS.NonHomRef.biallelic.VAF08.DP10.Norm.vcf.gz"

        list_path_gvcfs = [
            f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round3_SNV_Call/SNV_Calling/SNV_Calling.ONT.Hap{hap}.Triobin_{dict_hap_to_parental[hap]}/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.Polish_R3.SV.hap{hap}.ONT.triobin_{dict_hap_to_parental[hap]}.g.vcf.gz"
        ]

        path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round3_SNV_Call/SNV_Filtering/Triobinned_ONT_HomoAlt/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.Polish_R3.SV.hap{hap}.ONT.triobin_{dict_hap_to_parental[hap]}0a.PASS.NonHomRef.biallelic.VAF08.DP10.Norm.NAonSpecificGVCF.vcf"
        path_save_stat = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round3_SNV_Call/SNV_Filtering/Triobinned_ONT_HomoAlt/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.Polish_R3.SV.hap{hap}.ONT.triobin_{dict_hap_to_parental[hap]}0a.PASS.NonHomRef.biallelic.VAF08.DP10.Norm.NAonSpecificGVCF.vcf.stat"

        main(
            path_vcf,
            list_path_gvcfs,
            path_save,
            path_save_stat,
            n_jobs = 30
        )
        print(hap)
#%%
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
        for tgs in ["ONT", "HiFi"]:
            path_vcf = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round3_SNV_Call/SNV_Filtering/Triobinned_for_Evaluation/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.Polish_R3.SV.hap{hap}.{tgs}.triobin_{dict_hap_to_parental[hap]}0a.PASS.NonHomRef.Multiallelic_Norm.VAF05.NonHomRef.vcf.gz"

            list_path_gvcfs = [
                f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round3_SNV_Call/SNV_Calling/SNV_Calling.{tgs}.Hap{hap}.Triobin_{dict_hap_to_parental[hap]}/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.Polish_R3.SV.hap{hap}.{tgs}.triobin_{dict_hap_to_parental[hap]}.g.vcf.gz"
            ]

            path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round3_SNV_Call/SNV_Filtering/Triobinned_for_Evaluation/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.Polish_R3.SV.hap{hap}.{tgs}.triobin_{dict_hap_to_parental[hap]}0a.PASS.NonHomRef.Multiallelic_Norm.VAF05.NonHomRef.NAonSpecificGVCF.vcf"
            path_save_stat = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round3_SNV_Call/SNV_Filtering/Triobinned_for_Evaluation/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.Polish_R3.SV.hap{hap}.{tgs}.triobin_{dict_hap_to_parental[hap]}0a.PASS.NonHomRef.Multiallelic_Norm.VAF05.NonHomRef.NAonSpecificGVCF.vcf.stat"

            main(
                path_vcf,
                list_path_gvcfs,
                path_save,
                path_save_stat,
                n_jobs = 30
            )
            print(hap, tgs)
# %%
