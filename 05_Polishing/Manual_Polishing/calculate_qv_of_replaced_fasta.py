#%%
import sys,os,time
from pathlib import Path

sys.path.append(str(Path(__file__).parents[2] / "common"))
from Merqury_modified.run_qv_calculation_only import run as run_qv_calc



dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis"

curr_time = int(time.time())

k = 31

n_round = 7
afterround = 4

path_child_meryl = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Input_Filtering/PacBio/Meryl_Binning/PacBio_HiFi_ccs_converted.filt.Filter_L1000_QV20.k{k}.meryl"

list_ref_assemblies = ['Hifiasm_Binned_HiFi', 'Hifiasm_Binned_ALT_HiFi', 'Hifiasm_Binned_HiFi_with0', 'Hifiasm_Binned_ALT_HiFi_with0', 'Flye_Binned_ONT', 'Flye_Binned_ONT.with0', 'Flye_Binned_HiFi_with0', 'Flye_Binned_HiFi']

for hap in [1, 2]:
    for asm_from in list_ref_assemblies:
        path_asm = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Manual_Polishing/Round{n_round}_afterR{afterround}/Anchored_Sequences/Filtered_Sequences/Hap{hap}/Anchored_Sequence.{asm_from}.KOREF_S1.Trio_Specific.UpperCase.Polish_R{n_round-1}.Man.hap{hap}.sequence.without_anchor.fasta"
        # path_asm = path_asm_format

        assert os.path.exists(path_asm), path_asm

        dir_asm = os.path.dirname(path_asm)
        filename_asm = os.path.basename(path_asm)

        outputname_prefix = filename_asm.replace(".fasta", f".k{k}")

        path_output_prefix = os.path.join(dir_asm, "Merqury_QV", outputname_prefix)

        run_qv_calc(
            path_asm,
            path_child_meryl,
            path_output_prefix,
            f"{dir_workspace}/Polishing/Manual_Polishing/Merqury_Evaluation/tmp",
            f"{dir_workspace}/Polishing/Manual_Polishing/Merqury_Evaluation/sh",
            f"{dir_workspace}/Polishing/Manual_Polishing/Merqury_Evaluation/log",
            f"Evaluate_QV_only.Round{n_round}_afterR{afterround}.KOREF_S1.{asm_from}.Hap{hap}.k{k}.{curr_time}",
            # f"Evaluate_QV_only.KOREF_S1.{asm_name}.k21.{curr_time}",
            12,
            [],
            []
        )