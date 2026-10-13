#%%
import json, os, sys
from pathlib import Path

import pandas as pd
import numpy as np

sys.path.append(str(Path(__file__).parents[1]))
from fasta_parser import read_fasta_as_dict, save_dict_fasta_into_file, save_without_specific_seqnames_from_fasta_dict


def run(path_fasta, path_reassigned_order, path_save_prefix, manual_join, pseudo_gapsize = 100, n_pseudo_telomere_repeat = 5, telomere_repeat = "CCCTAA"):
    dict_reassigned_order = read_json(path_reassigned_order)

    path_save_joined_genome = f"{path_save_prefix}.reordered_scaffolds.n_pseudo_telorepeat{n_pseudo_telomere_repeat}.asm.fasta"
    path_save_leftover_fasta = f"{path_save_prefix}.unplaced.fasta"
    generate_joined_draft_genome(dict_reassigned_order, path_fasta, path_save_joined_genome, path_save_leftover_fasta, manual_join, pseudo_gapsize, n_pseudo_telomere_repeat, telomere_repeat)

def read_json(path_json):
    with open(path_json, 'r') as fr:
        dict_val = json.load(fr)
    return dict_val

def generate_joined_draft_genome(dict_reassigned_order, path_fasta, path_save_joined_genome, path_save_leftover_fasta, manual_join, pseudo_gapsize, n_pseudo_telomere_repeat, telomere_repeat):
    dict_fasta_input = read_fasta_as_dict(path_fasta)

    dict_draft_genome = dict()


    list_names_placed = list()
    for chrname, list_scaffold_order in dict_reassigned_order.items():
        manual_join_chr = manual_join.get(chrname, {})
        manual_join_chr_trim = manual_join_chr.get("trim", {})
        
        pseudo_gap = 'N' * manual_join_chr.get("gap", pseudo_gapsize)
        
        list_sequence = list()
        for ind_order, (scfname, ori) in enumerate(list_scaffold_order):
            if scfname == None:
                if ind_order == 0:
                    list_sequence.append(telomere_repeat * n_pseudo_telomere_repeat)
                elif ind_order == (len(list_scaffold_order)-1):
                    list_sequence.append(get_reverse_complement_sequence(telomere_repeat) * n_pseudo_telomere_repeat)
                else:
                    raise Exception(dict_reassigned_order)
            else:
                scf_sequence = dict_fasta_input[scfname].upper()
                assert ori in {'+', '-'}, dict_reassigned_order
                if ori == '-':
                    scf_sequence = get_reverse_complement_sequence(scf_sequence)
                    
                if manual_join_chr_trim.get(scfname) != None:
                    trim_left = manual_join_chr_trim[scfname]["trim_left"]
                    trim_right = manual_join_chr_trim[scfname]["trim_right"]
                    scf_sequence = scf_sequence[trim_left:len(scf_sequence)-trim_right]
                    
                list_sequence.append(scf_sequence)
                list_names_placed.append(scfname)
        dict_draft_genome[chrname] = pseudo_gap.join(list_sequence)
    
    save_dict_fasta_into_file(dict_draft_genome, path_save_joined_genome)
    save_without_specific_seqnames_from_fasta_dict(dict_fasta_input, list_names_placed, path_save_leftover_fasta)

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


#%%
if __name__ == "__main__":
    manual_join_per_hap = {
        1:{},
        2:{
            "chr14":{
                "trim":{
                    "h2tg000021l":{
                        "trim_left":0,
                        "trim_right":6971114
                    },
                    "ntLink_4.Left.0":{
                        "trim_left":0,
                        "trim_right":0
                    }
                },
                "gap":100                
            }
        }
    }
    
    
    for hap in [1, 2]:
        for n_pseudo_telomere_repeat in [0, 5, 10]:
            path_fasta_input = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling/Adjusted/Telomere_Detection/Fasta_Split_by_Telomere/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg.no_telo.telo_single_masked.fa.k32.w100.z1000.ntLink.scaffolds.adjusted_telomere.no_telo.telo_single.fasta"

            path_reassigned_order = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling.RagTag.default_option.CHM13/CuratedDraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.non_t2t_scaffold.hap{hap}.RagTag.default_option.minimap2_latest.CHM13/CuratedDraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.hap{hap}.RagTag_with_Telomere_Assign.minimap2_latest.CHM13.reordered_scaffolds.json"

            path_save_prefix = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling.RagTag.default_option.CHM13/CuratedDraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.non_t2t_scaffold.hap{hap}.RagTag.default_option.minimap2_latest.CHM13/CuratedDraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.hap{hap}.RagTag_with_Telomere_Assign.minimap2_latest.CHM13"

            run(path_fasta_input, path_reassigned_order, path_save_prefix, manual_join_per_hap[hap], pseudo_gapsize = 100, n_pseudo_telomere_repeat = n_pseudo_telomere_repeat, telomere_repeat = "CCCTAA")
# %%
