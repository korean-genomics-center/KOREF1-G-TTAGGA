#%%
import os, sys

import pandas as pd
from pathlib import Path

sys.path.append(str(Path(__file__).parents[1]))
from fasta_parser import read_fasta_as_dict, save_without_specific_seqnames_from_fasta_dict

#%%
hap = 2
refname = "CN1_combine"

dict_refname_to_path = {
    "CHM13" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/CHM13/chm13v2.0.fa",
    "CN1_combine" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/CN1/CN1_combine.v1.0.1.fasta"
}

dict_hap_to_sexchr_exclude = {
    1 : "chrX",
    2 : "chrY"
}

# list_path_tables_assigned = [
#     f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/Telomere_Detection/Fasta_Split_by_Telomere/Chromosome_Assign/Chromosome_Assign.KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.asm.dip.hap{hap}.p_ctg.telo_both.Ref_CHM13.tsv",
#     f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling/Adjusted/Telomere_Detection/Fasta_Split_by_Telomere/Chromosome_Assign/Chromosome_Assign.KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.asm.dip.hap{hap}.p_ctg.no_telo.telo_single_masked.fa.k32.w100.z1000.ntLink.scaffolds.adjusted_telomere.telo_both.Ref_CHM13.tsv"
# ]
# path_save_unassigned_ref = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling.RagTag.default_option.{refname}/{refname}.without_assigned_chromosomes_hap{hap}.fasta"


list_path_tables_assigned = [
    f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/Telomere_Detection/Fasta_Split_by_Telomere/Chromosome_Assign/Chromosome_Assign.KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg.telo_both.Ref_CHM13.tsv",
    f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling/Adjusted/Telomere_Detection/Fasta_Split_by_Telomere/Chromosome_Assign/Chromosome_Assign.KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg.no_telo.telo_single_masked.fa.k32.w100.z1000.ntLink.scaffolds.adjusted_telomere.telo_both.Ref_CHM13.tsv"
]
path_save_unassigned_ref = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling.RagTag.default_option.{refname}/{refname}.without_assigned_chromosomes_hap{hap}.fasta"

os.makedirs(os.path.dirname(path_save_unassigned_ref), exist_ok = True)
# %%
table_assigned = pd.concat(
    list(map(lambda path: pd.read_csv(path, sep = '\t'), list_path_tables_assigned))
)
# %%
def simplify_fasta_header(dict_fasta, sep = ' '):
    list_header = list(dict_fasta.keys())
    list_simplified_header = list(map(lambda val: val.split(sep)[0], list_header))

    assert len(list_simplified_header) == len(set(list_simplified_header)), list_simplified_header

    dict_header_to_simplified = dict(zip(list_header, list_simplified_header))

    dict_fasta_simplified = dict()
    for header, sequence in dict_fasta.items():
        dict_fasta_simplified[dict_header_to_simplified[header]] = sequence
    return dict_fasta_simplified

dict_ref_fasta = read_fasta_as_dict(dict_refname_to_path[refname])
dict_ref_fasta_simplified = simplify_fasta_header(dict_ref_fasta)

list_chrnames_assigned = table_assigned["Matched_Reference"].to_list()
chr_sex_exclude = dict_hap_to_sexchr_exclude[hap]
list_chrnames_assigned.append(chr_sex_exclude)

assert set(list_chrnames_assigned) < set(dict_ref_fasta_simplified.keys()), list_chrnames_assigned



save_without_specific_seqnames_from_fasta_dict(dict_ref_fasta_simplified, list_chrnames_assigned, path_save_unassigned_ref)

# %%
