#%%
import sys
from copy import deepcopy

sys.path.append("/BiO/Research/Project2/Endangered_Jindo_Dog/Resources/Scripts/KOGIC_KOREF_AGGA")
from fasta_parser import read_fasta_as_dict, save_dict_fasta_into_file

REVCOMP = {
    'A':'T',
    'T':'A',
    'G':'C',
    'C':'G'
}

def run(dict_fasta_paths, path_save):
    dict_filename_to_fasta = read_whole_fasta_files(dict_fasta_paths)
    dict_renamed_filename_to_fasta = rename_sequence_by_filename(dict_filename_to_fasta)
    
    dict_merged_fasta = dict()
    for dict_single_fasta in dict_renamed_filename_to_fasta.values():
        dict_merged_fasta = merge_another_fasta_file(dict_merged_fasta, dict_single_fasta)
    
    save_dict_fasta_into_file(dict_merged_fasta, path_save)

def read_whole_fasta_files(dict_fasta_paths):
    dict_filename_to_fasta = dict()
    for name, path_fasta in dict_fasta_paths.items():
        dict_filename_to_fasta[name] = read_fasta_as_dict(path_fasta)
    return dict_filename_to_fasta

def get_reverse_complement(sequence):
    assert set(sequence) == set(REVCOMP.keys()), set(sequence)
    
    rev_sequence = sequence[::-1]
    revcomp_sequence = ''.join(list(map(REVCOMP.__getitem__, rev_sequence)))
    return revcomp_sequence

def rename_sequence_by_filename(dict_filename_to_fasta):
    dict_renamed_filename_to_fasta = dict()
    for filename, dict_fasta in dict_filename_to_fasta.items():
        dict_renamed_filename_to_fasta[filename] = dict()
        for seqname, sequence in dict_fasta.items():
            dict_renamed_filename_to_fasta[filename][f"{filename}_{seqname}"] = sequence
    return dict_renamed_filename_to_fasta

def merge_another_fasta_file(dict_fasta, dict_additional_fasta):
    dict_merged_fasta = deepcopy(dict_fasta)
    
    list_sequences_exist = list()
    for sequence in dict_fasta.values():
        sequence_twice = sequence + sequence
        list_sequences_exist.append(sequence_twice)
        list_sequences_exist.append(get_reverse_complement(sequence_twice))
    
    for seqname, sequence in dict_additional_fasta.items():
        already_exists = False
        for sequence_exist in list_sequences_exist:
            if sequence in sequence_exist:
                already_exists = True
                break
        if not already_exists:
            dict_merged_fasta[seqname] = sequence
    return dict_merged_fasta


if __name__ == "__main2__":
    dict_fasta_paths = {}
    for k in [101, 151]:
        for ci in [10, 20]:
            filename = f"k{k}_ci{ci}"
            path_file = f"/BiO/Research/Project2/Endangered_Jindo_Dog/Results/Repeat_Annotation/Satellite_Repeat/run_srf/Hap2/Satellite_repeat_detection.J495799_Child.Hap2.k{k}.ci{ci}.srf_output.min_length_100.fa"
            dict_fasta_paths[filename] = path_file

    path_save = "/BiO/Research/Project2/Endangered_Jindo_Dog/Results/Repeat_Annotation/Satellite_Repeat/run_srf/Hap2/Satellite_repeat_detection.J495799_Child.Hap2.k101n151.ci10n20.srf_output.min_length_100.fa"
    run(dict_fasta_paths, path_save)
    

if __name__ == "__main__":
    dict_fasta_paths = {
        "Hap1":"/BiO/Research/Project2/Endangered_Jindo_Dog/Results/Repeat_Annotation/Satellite_Repeat/run_srf/Hap1/Satellite_repeat_detection.J495799_Child.Hap1.k101n151.ci10n20.srf_output.min_length_100.fa",
        "Hap2":"/BiO/Research/Project2/Endangered_Jindo_Dog/Results/Repeat_Annotation/Satellite_Repeat/run_srf/Hap2/Satellite_repeat_detection.J495799_Child.Hap2.k101n151.ci10n20.srf_output.min_length_100.fa"
    }
    path_save = "/BiO/Research/Project2/Endangered_Jindo_Dog/Results/Repeat_Annotation/Satellite_Repeat/run_srf/Mixed/Satellite_repeat_detection.J495799_Child.Hap1nHap2.k101n151.ci10n20.srf_output.min_length_100.fa"
    run(dict_fasta_paths, path_save)