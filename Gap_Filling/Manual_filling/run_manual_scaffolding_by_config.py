#%%
import sys,os,json,re
from pathlib import Path
from copy import deepcopy

sys.path.append(str(Path(__file__).parents[2]))
from fasta_parser import read_fasta_as_dict, save_dict_fasta_into_file


def run(path_config, path_save):
    dict_config = read_json(path_config)
    dict_sequence = read_input_sequence(dict_config["path_sequence"])

    dict_assembly = deepcopy(dict_sequence["assembly"])
    dict_modified_assembly = organize_contigs_by_manual_config(dict_assembly, dict_sequence, dict_config["manual_scaffolding"])

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

def organize_contigs_by_manual_config(dict_assembly, dict_sequence, dict_scaffolding_config):
    for chrname in dict_scaffolding_config.keys():
        list_contigs_chromosome, list_gaps_chromosome = split_sequence_by_gap(dict_assembly[chrname])

        list_sequences_organized = list()
        for ind, sequence_info in enumerate(dict_scaffolding_config[chrname]):
            if sequence_info["source"] == "assembly":
                sequence_append = list_contigs_chromosome[sequence_info["name_or_ind"]]
            else:
                sequence_append = dict_sequence[sequence_info["source"]][sequence_info["name_or_ind"]]
            sequence_append_modified = get_modified_sequence_from_sequence_info(sequence_append, sequence_info["strand"], sequence_info.get("range"))
            
            if ind > 0:
                gapseq = 'N'*100 if sequence_info["source"] != "assembly" else list_gaps_chromosome[sequence_info["name_or_ind"]-1]
                assert set(gapseq) == set('N')
                list_sequences_organized.append(gapseq)

            list_sequences_organized.append(sequence_append_modified)

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
    triobin = 'm'

    path_config = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon/Round5.TGS_Gapcloser_Modified_V2_re.Unpl_Scf.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.Manual_Scaff/Hap{hap}/scaffolding_config.json"
    path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon/Round5.TGS_Gapcloser_Modified_V2_re.Unpl_Scf.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.Manual_Scaff/Hap{hap}/CuratedDraftASM.KOREF_S1.Trio_Specific.hap{hap}.GapFilling.Unplaced_Scf.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.Manual_Scaff.fasta"

    run(path_config, path_save)
# %%
