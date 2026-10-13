#%%
import sys, os
from pathlib import Path

import pandas as pd

sys.path.append(str(Path(__file__).parents[2] / "common"))
from fasta_parser import read_fasta_as_dict, save_dict_fasta_into_file


def rename_teloboth_fasta_as_chromosome(path_teloboth, path_chromosome_assign, path_save):
    os.makedirs(os.path.dirname(path_save), exist_ok = True)

    dict_fasta = read_fasta_as_dict(path_teloboth)
    table_assign = pd.read_csv(path_chromosome_assign, sep = '\t')

    dict_assigned_fasta = dict()
    for _, row in table_assign.iterrows():
        sequence = dict_fasta[row["Scaffold"]].upper()
        if row["Strand"] == '-':
            sequence = get_reverse_complement_sequence(sequence)
        dict_assigned_fasta[row["Matched_Reference"]] = sequence

    save_dict_fasta_into_file(dict_assigned_fasta, path_save)

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
    for hap in [1, 2]:
        path_teloboth = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/Telomere_Detection/Fasta_Split_by_Telomere/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg.telo_both.fasta"
        path_chromosome_assign = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/Telomere_Detection/Fasta_Split_by_Telomere/Chromosome_Assign/Chromosome_Assign.KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg.telo_both.Ref_CHM13.tsv"

        path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling.RagTag.default_option.CHM13/CuratedDraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.non_t2t_scaffold.hap{hap}.RagTag.default_option.minimap2_latest.CHM13/CuratedDraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.telo_both.hap{hap}.assigned_CHM13.asm.fasta"

        rename_teloboth_fasta_as_chromosome(path_teloboth, path_chromosome_assign, path_save)
# %%
if __name__ == "__main__":
    for hap in [1, 2]:
        path_teloboth = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling/Adjusted/Telomere_Detection/Fasta_Split_by_Telomere/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg.no_telo.telo_single_masked.fa.k32.w100.z1000.ntLink.scaffolds.adjusted_telomere.telo_both.fasta"
        path_chromosome_assign = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling/Adjusted/Telomere_Detection/Fasta_Split_by_Telomere/Chromosome_Assign/Chromosome_Assign.KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg.no_telo.telo_single_masked.fa.k32.w100.z1000.ntLink.scaffolds.adjusted_telomere.telo_both.Ref_CHM13.tsv"

        path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling.RagTag.default_option.CHM13/CuratedDraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.non_t2t_scaffold.hap{hap}.RagTag.default_option.minimap2_latest.CHM13/CuratedDraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.telo_both.hap{hap}.assigned_CHM13.asm.fasta"

        rename_teloboth_fasta_as_chromosome(path_teloboth, path_chromosome_assign, path_save)
# %%
