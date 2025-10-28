#%%
import sys,os,time
from pathlib import Path

sys.path.append(str(Path(__file__).parents[2]))
from Merqury_modified.run_hapmer_count import run as run_hamming_calc



dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis"

curr_time = int(time.time())

k = 31

n_round = 5
afterround = 2

path_hapmer_paternal = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/PacBio/Hapmer_Diff/Hapmer_diff_k{k}/PATERNAL_Trimmed_F15_T5.noNbase.merged.k{k}.hapmer.meryl"
path_hapmer_maternal = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/PacBio/Hapmer_Diff/Hapmer_diff_k{k}/MATERNAL_Trimmed_F15_T5.noNbase.merged.k{k}.hapmer.meryl"

list_ref_assemblies = ['Hifiasm_Binned_HiFi', 'Hifiasm_Binned_ALT_HiFi', 'Hifiasm_Binned_HiFi_with0', 'Hifiasm_Binned_ALT_HiFi_with0', 'Flye_Binned_ONT', 'Flye_Binned_ONT.with0', 'Flye_Binned_HiFi_with0', 'Flye_Binned_HiFi']

for hap in [1, 2]:
    for asm_from in list_ref_assemblies:
        path_asm = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Manual_Polishing/Round{n_round}_afterR{afterround}/Anchored_Sequences/Filtered_Sequences/Hap{hap}/Anchored_Sequence.{asm_from}.KOREF_S1.Trio_Specific.UpperCase.Polish_R{n_round-1}.Man.hap{hap}.sequence.with_anchor.fasta"
        # path_asm = path_asm_format

        assert os.path.exists(path_asm), path_asm

        dir_asm = os.path.dirname(path_asm)
        filename_asm = os.path.basename(path_asm)

        outputname_prefix = filename_asm.replace(".fasta", '')

        path_output_prefix = os.path.join(dir_asm, "Merqury_Hamming", outputname_prefix)

        run_hamming_calc(
            path_asm,
            path_hapmer_paternal,
            path_hapmer_maternal,
            path_output_prefix,
            f"{dir_workspace}/Polishing/Manual_Polishing/Merqury_Evaluation/tmp",
            f"{dir_workspace}/Polishing/Manual_Polishing/Merqury_Evaluation/sh",
            f"{dir_workspace}/Polishing/Manual_Polishing/Merqury_Evaluation/log",
            f"Evaluate_Hamming_only.Round{n_round}_afterR{afterround}.KOREF_S1.{asm_from}.Hap{hap}.k{k}.{curr_time}",
            32,
            [],
            []
        )