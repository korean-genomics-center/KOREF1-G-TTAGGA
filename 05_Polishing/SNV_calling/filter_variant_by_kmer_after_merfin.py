#%%
from pathlib import Path
import os,sys
from itertools import combinations, chain

from joblib import Parallel, delayed
import pandas as pd
import numpy as np
from Bio.bgzf import BgzfReader

sys.path.append(str(Path(__file__).parents[2] / "common"))
from fasta_parser import read_fasta_as_dict


def run(path_vcf, path_reference, path_kmer_count, path_output, kmer_size, n_jobs):
    table_vcf, vcf_header = read_vcf_as_table(path_vcf, return_header=True)
    dict_ref = read_fasta_as_dict(path_reference)
    
    dict_grouped_variants = group_variants_by_kmer_size(table_vcf, kmer_size)
    list_all_kmers_check = list(set(get_all_kmers_for_check_from_variant_combination(dict_grouped_variants, kmer_size, dict_ref)))
    list_all_kmers_check_rc = list(map(get_reverse_complementary, list_all_kmers_check))
    
    set_kmers_check = set(list_all_kmers_check + list_all_kmers_check_rc)
    
    set_readmers = read_kmer_count_parallel(path_kmer_count, set_kmers_check, n_jobs)
    readmers_rv = list(map(get_reverse_complementary, set_readmers))
    
    all_readmers = set(list(set_readmers)+readmers_rv)
    
    table_vcf_best_variants = get_best_variant_combinations_to_reduce_error(dict_grouped_variants, all_readmers, kmer_size, dict_ref)
    save_vcf_with_header(table_vcf_best_variants, vcf_header, path_output)

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
    
def read_kmer_count(path_kmer_count):
    list_kmers = list()
    
    with open(path_kmer_count, 'r') as fr:
        for line in fr:
            sequence, cnt = line.strip('\n').split()
            list_kmers.append(sequence)
    
    set_kmers = set(list_kmers)
    return set_kmers
    
def read_kmer_count_parallel(path_kmer_count, set_kmers_check, n_jobs = 20):
    filesize = os.path.getsize(path_kmer_count)
    list_seekpoints = np.linspace(0, filesize, n_jobs + 1)
    
    list_seekpoints_readline = [list_seekpoints[0]]
    with open(path_kmer_count, 'r') as fr:
        for seekpoint in list_seekpoints[1:-1]:
            fr.seek(seekpoint)
            fr.readline()
            list_seekpoints_readline.append(fr.tell())
    list_seekpoints_readline.append(filesize)
    
    def read_kmer_count_in_range(seek_start, seek_end):
        list_kmers_part = list()
        file_reader = open(path_kmer_count, 'r')
        file_reader.seek(seek_start)
        
        while 1:
            line = file_reader.readline()
            if not line:
                break
            sequence, cnt = line.strip('\n').split()
            if sequence in set_kmers_check:
                list_kmers_part.append(sequence)
            if file_reader.tell() >= seek_end:
                break
        
        file_reader.close()
        return list_kmers_part
    
    with Parallel(n_jobs = n_jobs) as parallel:
        list_of_list_kmers = parallel(
            delayed(read_kmer_count_in_range)(
                list_seekpoints_readline[ind],
                list_seekpoints_readline[ind+1]
            ) for ind in range(len(list_seekpoints_readline)-1)
        )
    
    list_kmers = list(chain(*list_of_list_kmers))
    
    set_kmers = set(list_kmers)
    return set_kmers

def group_variants_by_kmer_size(table_vcf, kmer_size):
    dict_chrname_to_variant_group = dict()
    
    for chrname in table_vcf["#CHROM"].unique():
        dict_chrname_to_variant_group[chrname] = list()
        
        table_vcf_chrname = table_vcf[table_vcf["#CHROM"] == chrname]
        
        list_pos_on_vcf = table_vcf_chrname["POS"].to_list()
        
        list_groups = list()
        group = list()
        for pos in sorted(list_pos_on_vcf):
            if len(group) == 0:
                group.append(pos)
            if pos - group[-1] < kmer_size:
                group.append(pos)
            else:
                list_groups.append(group)
                group = [pos]
        if len(group) > 0:
            list_groups.append(group)
        
        for group in list_groups:
            dict_chrname_to_variant_group[chrname].append(table_vcf_chrname[table_vcf_chrname["POS"].isin(group)])
    return dict_chrname_to_variant_group

def get_all_kmers_for_check_from_variant_combination(dict_grouped_variants, kmer_size, dict_fasta):
    list_all_kmers_check = list()
    
    for chrname, list_tables_by_group in dict_grouped_variants.items():
        for table_vcf in list_tables_by_group:
            list_kmers_from_vcf = collect_kmers_for_check_from_vcf(table_vcf, kmer_size, dict_fasta[chrname])
            list_all_kmers_check.extend(list_kmers_from_vcf)
    return list_all_kmers_check

def collect_kmers_for_check_from_vcf(table_vcf, kmer_size, ref_sequence):
    min_pos = table_vcf["POS"].min()
    max_pos = max(table_vcf["POS"] + table_vcf["REF"].apply(len) - 1)
    
    pos_collect_start = max(1, min_pos - kmer_size + 1)
    pos_collect_end = min(max_pos + kmer_size - 1, len(ref_sequence))
    
    table_vcf_sorted = table_vcf.sort_values(by = ["POS"]).reset_index(drop = True)
    sequence_collect = ref_sequence[pos_collect_start-1:pos_collect_end]
    
    list_kmers_check = get_canonical_kmer_from_sequence(sequence_collect, kmer_size)
    
    list_variant_indexes = list(table_vcf_sorted.index)
    list_variant_combinations = get_all_combinations(list_variant_indexes)
    
    for variant_comb in list_variant_combinations:
        sequence_from_variant_comb = get_variant_replaced_sequence(table_vcf_sorted.loc[variant_comb, :], pos_collect_start, sequence_collect)
        list_kmers_check.extend(get_canonical_kmer_from_sequence(sequence_from_variant_comb, kmer_size))
    return list_kmers_check        

def get_best_variant_combinations_to_reduce_error(dict_grouped_variants, set_readmers, kmer_size, dict_fasta):
    list_tables_best_comb = list()
    
    for chrname, list_tables_by_group in dict_grouped_variants.items():
        for table_vcf in list_tables_by_group:
            table_vcf_best_comb = check_best_combinations_of_variant_group(table_vcf, set_readmers, kmer_size, dict_fasta[chrname])
            list_tables_best_comb.append(table_vcf_best_comb)
    return pd.concat(list_tables_best_comb)


def check_best_combinations_of_variant_group(table_vcf, set_readmers, kmer_size, ref_sequence):
    min_pos = table_vcf["POS"].min()
    max_pos = max(table_vcf["POS"] + table_vcf["REF"].apply(len) - 1)
    
    pos_collect_start = max(1, min_pos - kmer_size + 1)
    pos_collect_end = min(max_pos + kmer_size - 1, len(ref_sequence))
    
    table_vcf_sorted = table_vcf.sort_values(by = ["POS"]).reset_index(drop = True)
    sequence_collect = ref_sequence[pos_collect_start-1:pos_collect_end]
    
    list_variant_indexes = list(table_vcf_sorted.index)
    list_variant_combinations = get_all_combinations(list_variant_indexes)
    
    n_novelkmer_original = count_unobserved_kmer_count_of_sequence(sequence_collect, set_readmers, kmer_size)
    
    dict_comb_to_sequence = dict()
    dict_comb_to_novelkmer = dict()
    for variant_comb in list_variant_combinations:
        sequence_from_variant_comb = get_variant_replaced_sequence(table_vcf_sorted.loc[variant_comb, :], pos_collect_start, sequence_collect)
        dict_comb_to_sequence[variant_comb] = sequence_from_variant_comb
        dict_comb_to_novelkmer[variant_comb] = count_unobserved_kmer_count_of_sequence(sequence_from_variant_comb, set_readmers, kmer_size)
    
    list_comb_to_novelkmer = [((), n_novelkmer_original)] + list(dict_comb_to_novelkmer.items())
    
    min_novelkmer = sorted(list_comb_to_novelkmer, key = lambda val: val[1])[0][1]
    
    list_comb_with_minimum_novelkmer = list(filter(lambda val: val[1] == min_novelkmer, list_comb_to_novelkmer))
    
    if len(list_comb_with_minimum_novelkmer) == 1 and len(list_comb_with_minimum_novelkmer[0][0]) == 0:
        return pd.DataFrame()
    
    maxvariants_with_minimum_novelkmer = sorted(list_comb_with_minimum_novelkmer, key = lambda val: len(val[0]), reverse = True)[0]
    
    table_vcf_of_best = table_vcf_sorted.loc[maxvariants_with_minimum_novelkmer[0], :]
    return table_vcf_of_best    
    
def get_all_combinations(list_values):
    list_comb = list()
    
    for i in range(1, len(list_values)+1):
        list_comb.extend(list(combinations(list_values, i)))
    return list_comb    

def get_variant_replaced_sequence(table_vcf, pos_collect_start, sequence_collect):
    sequence_modified = sequence_collect
    
    for _, row_variant in table_vcf.sort_values(by = ["POS"], ascending = False).iterrows():
        pos_variant = row_variant["POS"]
        ref_variant = row_variant["REF"]
        alt_variant = row_variant["ALT"]
        
        ind_variant_in_collected = pos_variant - pos_collect_start
        
        sequence_modified_front = sequence_modified[:ind_variant_in_collected]
        sequence_modified_end = sequence_modified[ind_variant_in_collected+len(ref_variant):]
        sequence_original = sequence_collect[ind_variant_in_collected:ind_variant_in_collected+len(ref_variant)]
        assert sequence_original == ref_variant, f"{pos_variant} : {sequence_original} > {ref_variant}"
        
        sequence_modified = f"{sequence_modified_front}{alt_variant}{sequence_modified_end}"
    return sequence_modified

def get_canonical_kmer_from_sequence(sequence, kmer_size):
    list_kmers = list()
    
    for ind in range(len(sequence) - kmer_size):
        kmer_sequence = sequence[ind:ind+kmer_size]
        kmer_sequence_revcomp = get_reverse_complementary(kmer_sequence)
        
        if kmer_sequence < kmer_sequence_revcomp:
            list_kmers.append(kmer_sequence)
        else:
            list_kmers.append(kmer_sequence_revcomp)
    return list_kmers

def get_reverse_complementary(sequence):
    dict_rev = {
        'A':'T',
        'T':'A',
        'G':'C',
        'C':'G'
    }
    
    return ''.join(list(map(dict_rev.__getitem__, sequence[::-1])))

def count_unobserved_kmer_count_of_sequence(sequence, set_readmer, kmer_size):
    list_kmer_of_sequence = get_canonical_kmer_from_sequence(sequence, kmer_size)
    list_nonexistence_of_kmer = list(map(lambda val: val not in set_readmer, list_kmer_of_sequence))
    
    return sum(list_nonexistence_of_kmer)

def save_vcf_with_header(table_vcf, header, path_output):
    with open(path_output, 'w') as fw:
        fw.write(header)
        fw.write(get_context_to_write(table_vcf.columns))
        for _, row in table_vcf.iterrows():
            fw.write(get_context_to_write(row.to_list()))

def get_context_to_write(list_values, sep = '\t'):
    return sep.join(list(map(str, list_values))) + '\n'
    


# %%
if __name__ == "__main__":
    
    for hap in [1, 2]:
        kmer_size = 31
        
        path_vcf = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round8_SNV_Call/Polishing_Target/Merfin_Filtered/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.Polish_R8.SV.hap{hap}.HiFi.triobin_pm0a.ONT_triobin_specific_HomAlt_VAF08DP10.PASS.NonHomRef.ForPolishing.Merged_AllTarget.sorted.polish.vcf"
        path_reference = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Reference_Fasta/Round8.SV_Polished_Fasta/PolishedASM.KOREF_S1.Trio_Specific.UpperCase.Polish_R8.SV.hap{hap}.fasta"
        path_kmer_count = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/PacBio/Meryl_Binning/PacBio_HiFi_ccs_converted.filt.Filter_L1000_QV20.k{kmer_size}.gt1.kmer_count.txt"
        
        
        n_jobs = 90
        
        path_output = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round8_SNV_Call/Polishing_Target/Merfin_Filtered/SNV_Call.PepperDeepvariant.KOREF_PolishedASM.Polish_R8.SV.hap{hap}.HiFi.triobin_pm0a.ONT_triobin_specific_HomAlt_VAF08DP10.PASS.NonHomRef.ForPolishing.Merged_AllTarget.sorted.polish.best_k{kmer_size}.vcf"
        
        # run(path_vcf, path_reference, path_kmer_count, path_output, kmer_size, n_jobs)
# %%
