#%%
import sys,os,json,re
from pathlib import Path
from copy import deepcopy

sys.path.append("/BiO/Research/Project2/Endangered_Goral/Resources/Scripts/KOGIC_KOREF_AGGA")
from fasta_parser import read_fasta_as_dict, save_dict_fasta_into_file


def run(path_config, path_save):
    dict_config = read_json(path_config)
    dict_sequence = read_input_sequence(dict_config["path_sequence"])

    dict_modified_assembly = organize_contigs_by_manual_config(dict_sequence, dict_config["manual_scaffolding"])

    save_dict_fasta_into_file(dict_modified_assembly, path_save)

def read_json(path_json):
    with open(path_json, 'r') as fr:
        dict_val = json.load(fr)
    return dict_val

def read_input_sequence(dict_path_sequence):
    dict_sequence = dict()
    dict_sequence["assembly"] = read_fasta_as_dict(dict_path_sequence["assembly"])
    for seqname, path_sequence in dict_path_sequence["scaffolding_contig"].items():
        assert dict_sequence.get(seqname) == None, f"Name of input contig overlaps: {seqname}"
        dict_sequence[seqname] = read_fasta_as_dict(path_sequence)
    return dict_sequence

def organize_contigs_by_manual_config(dict_sequence, dict_scaffolding_config):
    dict_assembly = dict()
    
    for chrname in dict_scaffolding_config.keys():
        print(chrname)
        list_sequences_organized = list()
        for ind, sequence_info in enumerate(dict_scaffolding_config[chrname]):
            sequence_append = dict_sequence[sequence_info["source"]][sequence_info["name_or_ind"]]
            sequence_append_modified = get_modified_sequence_from_sequence_info(sequence_append, sequence_info["strand"], sequence_info.get("range"))
            
            if ind > 0:
                gapseq = 'N'*100
                assert set(gapseq) == set('N')
                list_sequences_organized.append(gapseq)

            list_sequences_organized.append(sequence_append_modified)
        
        from collections import Counter
        for sequence in list_sequences_organized:
            print(Counter(sequence))

        dict_assembly[chrname] = ''.join(list_sequences_organized)
    return dict_assembly

def split_sequence_by_gap(sequence):
    re_findgap = re.compile("[N]+")

    sequence_upp = sequence.upper()
    list_contigs = list(filter(lambda val: len(val) > 0, sequence_upp.split('N')))

    list_gaps = list()
    for gap_found in re_findgap.finditer(sequence_upp):
        gap_region = gap_found.span()
        list_gaps.append(sequence_upp[gap_region[0]:gap_region[1]])
    return (list_contigs, list_gaps)

def get_modified_sequence_from_sequence_info(sequence, strand, range):
    if range != None:
        sequence_in_range = sequence[range[0]-1:range[1]]
    else:
        sequence_in_range = sequence
    if strand == '-':
        sequence_in_range = get_reverse_complement_sequence(sequence_in_range)
    return sequence_in_range
    
def get_reverse_complement_sequence(sequence):
    dict_reverse = {
        'A' : 'T',
        'T' : 'A',
        'G' : 'C',
        'C' : 'G',
        'N' : 'N'
    }

    reversed_sequence = sequence[::-1]
    reverse_complement_sequence = list(map(dict_reverse.__getitem__, reversed_sequence))
    return ''.join(reverse_complement_sequence)

# %%
if __name__ == "__main__":
    hap = 2

    path_config = f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Manual_Scaffold/scaffolding_config.json"
    path_save = f"/BiO/Research/Project2/Endangered_Goral/Results/Genome_Assembly/Assembly_Scaffolding/Hifiasm/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.default_option/Manual_Scaffold/KoGoral.ont_L50kb_QV15.hifi_1kb_QV20.with_HiC.asm.hic.manual_scaffold.fasta"

    run(path_config, path_save)
# %%
