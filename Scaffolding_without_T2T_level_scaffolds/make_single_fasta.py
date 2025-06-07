#%%
from pathlib import Path
import sys,os

sys.path.append(str(Path(__file__).parents[1]))
from fasta_parser import *


#%%
############### Reference
list_save_ref = ["chr13", "chr15"]

dict_fasta = read_fasta_as_dict("/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/CN1/CN1_combine.v1.0.1.fasta")
# dict_fasta = read_fasta_as_dict("/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/CHM13/chm13v2.0.fa")
#%%
prefix_save = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling.RagTag.default_option.CHM13/DraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.telo_single.hap2.RagTag.default_option.minimap2_latest.CHM13/OneToOneMapping/Raw_Fasta/Reference/CN1_combine"

for chrname, sequence in dict_fasta.items():
    chrname = chrname.split()[0]

    if chrname in list_save_ref:
        len_seq = len(sequence)
        # middle_sequence = int(len_seq/2)
        middle_sequence = int(1e6)
        sequence_front = sequence[:middle_sequence]
        sequence_end = sequence[-middle_sequence:]


        dict_seq = {f"{chrname}":sequence}
        dict_seq_front = {f"{chrname}.Front1Mb":sequence_front}
        dict_seq_end = {f"{chrname}.End1Mb":sequence_end}

        save_dict_fasta_into_file(
            dict_seq,
            f"{prefix_save}.{chrname}.fasta"
        )
        save_dict_fasta_into_file(
            dict_seq_front,
            f"{prefix_save}.{chrname}.Front1Mb.fasta"
        )
        save_dict_fasta_into_file(
            dict_seq_end,
            f"{prefix_save}.{chrname}.End1Mb.fasta"
        )

#%%
####################### QRY
hap = 1

list_save_qry = ["h1tg000060l", "h1tg000044l", "h1tg000041l", "h1tg000032l"]

dict_fasta = read_fasta_as_dict(f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling/Adjusted/Telomere_Detection/Fasta_Split_by_Telomere/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap{hap}.p_ctg.no_telo.telo_single_masked.fa.k32.w100.z1000.ntLink.scaffolds.adjusted_telomere.telo_single.fasta")

dir_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/ntLink.telomasked.without_T2T_scaffold.ONT_triobin_specific.sensitive_gapfilling.RagTag.default_option.CHM13/DraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.telo_single.hap{hap}.RagTag.default_option.minimap2_latest.CHM13/OneToOneMapping/Raw_Fasta/Scaffolds"

for scfname, sequence in dict_fasta.items():
    if scfname in list_save_qry:
        dict_seq = {f"{scfname}":sequence}

        save_dict_fasta_into_file(
            dict_seq,
            f"{dir_save}/{scfname}.fasta"
        )
# %%
