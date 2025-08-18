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
    list_path_gvcf_difftrio,
    path_save
):
    table_vcf, vcf_header = read_vcf_as_table(path_vcf, True)
    list_vcfs_sametrio = list(map(read_vcf_as_table, list_path_vcf_sametrio))
    list_vcfs_difftrio = list(map(read_vcf_as_table, list_path_vcf_difftrio))
    list_gvcfs_region_difftrio = list(map(read_gvcf_region, list_path_gvcf_difftrio))

    table_vcf["VarID"] = get_varids_from_vcf(table_vcf)

    dict_varid_to_gtset_sametrio = get_snvs_with_genotype_combinations_from_vcfs(list_vcfs_sametrio)
    dict_varid_to_gtset_difftrio = get_snvs_with_genotype_combinations_from_vcfs(list_vcfs_difftrio)

    list_varid_to_polish, list_varid_to_ambiguous = extract_polishing_target_snv_from_sametrio_call(table_vcf["VarID"], dict_varid_to_gtset_sametrio)

    list_varid_to_ambiguous_notcalled_in_difftrio = list(set(list_varid_to_ambiguous) - set(dict_varid_to_gtset_difftrio.keys()))
    list_varid_amb_evaluated_on_gvcf = evaluate_ambiguous_varids_on_difftrio_gvcf(list_varid_to_ambiguous_notcalled_in_difftrio, list_gvcfs_region_difftrio)

    list_total_varids_for_polishing = list_varid_to_polish+list_varid_amb_evaluated_on_gvcf

    table_vcf_for_polishing = table_vcf[table_vcf["VarID"].isin(list_total_varids_for_polishing)]
    table_vcf_for_polishing = table_vcf_for_polishing.drop(columns = ["VarID"])

    save_vcf_for_polishing(table_vcf_for_polishing, vcf_header, path_save)


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

def extract_polishing_target_snv_from_sametrio_call(list_varids_input, dict_varid_to_gtset_sametrio):
    dict_sametrio_gtset_to_input_vcf_varids, list_varids_notin_gtset = classify_varids_by_genotype_combination(list_varids_input, dict_varid_to_gtset_sametrio)

    list_gtsets_with_homalt, list_gtsets_without_homalt = classify_tuple_whether_containing_specific_value(dict_sametrio_gtset_to_input_vcf_varids.keys(), "1/1")

    list_direct_polishing = list(chain(*list(map(dict_sametrio_gtset_to_input_vcf_varids.__getitem__, list_gtsets_with_homalt))))
    list_ambiguous_polishing = list(chain(*list(map(dict_sametrio_gtset_to_input_vcf_varids.__getitem__, list_gtsets_without_homalt))))

    return list_direct_polishing, list_ambiguous_polishing

def classify_varids_by_genotype_combination(list_varid, dict_varid_to_gtset):
    dict_gtset_to_varids = dict()

    for varid in list_varid:
        gtset = dict_varid_to_gtset.get(varid, "None")
        if dict_gtset_to_varids.get(gtset) == None:
            dict_gtset_to_varids[gtset] = list()
        dict_gtset_to_varids[gtset].append(varid)
    list_varids_nonexist = dict_gtset_to_varids.pop("None")
    return dict_gtset_to_varids, list_varids_nonexist

def classify_tuple_whether_containing_specific_value(list_tuples, value_detect):
    list_tuples_with_value = list(filter(lambda tup: value_detect in tup, list_tuples))
    list_tuples_without_value = list(filter(lambda tup: value_detect not in tup, list_tuples))
    return list_tuples_with_value, list_tuples_without_value

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

def evaluate_ambiguous_varids_on_difftrio_gvcf(list_varids_check, list_gvcfs_region_difftrio):
    dict_varids_to_gts_on_gvcf = check_filtered_ambiguous_snvs_variant_calls_on_gvcf(list_varids_check, list_gvcfs_region_difftrio)

    list_varids_homoref_on_gvcf = list()
    for varid, list_gts_on_gvcf in dict_varids_to_gts_on_gvcf.items():
        set_gts_on_gvcf = set(list_gts_on_gvcf)
        if {"0/0"} <= set_gts_on_gvcf and set_gts_on_gvcf <= {"0/0", "./."}:
            list_varids_homoref_on_gvcf.append(varid)
    return list_varids_homoref_on_gvcf

def check_filtered_ambiguous_snvs_variant_calls_on_gvcf(list_varids_check, list_gvcfs_region_difftrio):
    for dict_gvcf_region in list_gvcfs_region_difftrio:
        verify_gvcf_regions_sorted(dict_gvcf_region)
    
    dict_varid_to_gvcfcall = dict()
    for varid in list_varids_check:
        dict_varid_to_gvcfcall[varid] = list()
        varid_chrname, varid_pos, _, _ = varid.split(':')
        varid_pos = int(varid_pos)

        for ind, dict_gvcf_region in enumerate(list_gvcfs_region_difftrio):
            list_regions_on_gvcf = list(dict_gvcf_region[varid_chrname].keys())
            ind_region_containing_position = search_putative_index_region_containing_position(
                list_regions_on_gvcf,
                varid_pos, 0, len(list_regions_on_gvcf)-1
            )
            assert ind_region_containing_position != None
            list_ranges_ind_check = range(max(0, ind_region_containing_position-10), min(ind_region_containing_position+11, len(list_regions_on_gvcf)-1))
            list_regions_check = list(map(lambda ind: list_regions_on_gvcf[ind], list_ranges_ind_check))
            list_regions_containing_pos = list(filter(lambda region: check_is_value_in_range(region, varid_pos), list_regions_check))
            assert list_regions_check[0] != list_regions_containing_pos[0] and list_regions_check[-1] != list_regions_containing_pos[-1]

            list_gts_containing_pos = list(map(dict_gvcf_region[varid_chrname].__getitem__, list_regions_containing_pos))

            assert len(list_gts_containing_pos) > 0

            dict_varid_to_gvcfcall[varid].extend(list_gts_containing_pos)
    return dict_varid_to_gvcfcall

def verify_gvcf_regions_sorted(dict_gvcf_region):
    for chrname in dict_gvcf_region.keys():
        list_regions = list(dict_gvcf_region[chrname])

        region_before = list_regions[0]
        for region in list_regions[1:]:
            assert region[0] > region_before[0], f"{chrname}, {region_before}, {region}"
            region_before = region

def search_putative_index_region_containing_position(list_regions, pos_find, ind_searchblock_start, ind_searchblock_end):
    ind_searchblock_middle = int((ind_searchblock_start+ind_searchblock_end)/2)
    region_middle = list_regions[ind_searchblock_middle]
    region_start = list_regions[ind_searchblock_start]
    region_end = list_regions[ind_searchblock_end]

    if ind_searchblock_start >= ind_searchblock_end:
        if check_is_value_in_range(region_middle, pos_find):
            return region_middle
        else:
            return None

    if check_is_value_in_range(region_middle, pos_find):
        return ind_searchblock_middle
    elif check_is_value_in_range(region_start, pos_find):
        return ind_searchblock_start
    elif check_is_value_in_range(region_end, pos_find):
        return ind_searchblock_end
    elif region_middle[0] > pos_find:
        return search_putative_index_region_containing_position(list_regions, pos_find, ind_searchblock_start, ind_searchblock_middle-1)
    else:
        return search_putative_index_region_containing_position(list_regions, pos_find, ind_searchblock_middle+1, ind_searchblock_end)

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
    dict_hap_to_triobin = {
        1:'p',
        2:'m'
    }
    dict_triobin_exchange = {
        'p':'m',
        'm':'p'
    }

    for hap in [2]:
        triobin_same = dict_hap_to_triobin[hap]
        triobin_diff = dict_triobin_exchange[triobin_same]

        path_vcf = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Align_Polished_ASM.From_R1.SV/SNV_Call/SNV_Calling.HiFi.Hap{hap}.Triobin_pm0a/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.R1_Sv.hap{hap}.HiFi.triobin_pm0a.PASS.NonHomRef.biallelic.GQ30.VAF04.Norm.Het.vcf.gz"
        path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Align_Polished_ASM.From_R1.SV/SNV_Call/SNV_Calling.HiFi.Hap{hap}.Triobin_pm0a/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.R1_Sv.hap{hap}.HiFi.triobin_pm0a.PASS.NonHomRef.biallelic.GQ30.VAF04.Norm.Het.ForPolishing.vcf"

        list_path_vcf_sametrio = [
            f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Align_Polished_ASM.From_R1.SV/SNV_Call/SNV_Calling.HiFi.Hap{hap}.Triobin_{triobin_same}/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.R1_Sv.hap{hap}.HiFi.triobin_{triobin_same}.PASS.NonHomRef.Multiallelic_Norm.VAF05.NonHomRef.vcf.gz",
            f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Align_Polished_ASM.From_R1.SV/SNV_Call/SNV_Calling.ONT.Hap{hap}.Triobin_{triobin_same}/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.R1_Sv.hap{hap}.ONT.triobin_{triobin_same}.PASS.NonHomRef.Multiallelic_Norm.VAF05.NonHomRef.vcf.gz"
        ]
        list_path_vcf_difftrio = [
            f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Align_Polished_ASM.From_R1.SV/SNV_Call/SNV_Calling.HiFi.Hap{hap}.Triobin_{triobin_diff}/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.R1_Sv.hap{hap}.HiFi.triobin_{triobin_diff}.PASS.NonHomRef.Multiallelic_Norm.VAF05.NonHomRef.vcf.gz",
            f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Align_Polished_ASM.From_R1.SV/SNV_Call/SNV_Calling.ONT.Hap{hap}.Triobin_{triobin_diff}/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.R1_Sv.hap{hap}.ONT.triobin_{triobin_diff}.PASS.NonHomRef.Multiallelic_Norm.VAF05.NonHomRef.vcf.gz"
        ]
        list_path_gvcf_difftrio = [
            f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Align_Polished_ASM.From_R1.SV/SNV_Call/SNV_Calling.HiFi.Hap{hap}.Triobin_{triobin_diff}/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.R1_Sv.hap{hap}.HiFi.triobin_{triobin_diff}.g.vcf.gz",
            f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Align_Polished_ASM.From_R1.SV/SNV_Call/SNV_Calling.ONT.Hap{hap}.Triobin_{triobin_diff}/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.R1_Sv.hap{hap}.ONT.triobin_{triobin_diff}.g.vcf.gz"
        ]

        # main(
        #     path_vcf,
        #     list_path_vcf_sametrio,
        #     list_path_vcf_difftrio,
        #     list_path_gvcf_difftrio,
        #     path_save
        # )
# %%
