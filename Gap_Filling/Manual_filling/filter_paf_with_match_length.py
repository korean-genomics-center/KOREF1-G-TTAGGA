#%%
import os,sys
import subprocess

from joblib import Parallel, delayed

def run(
    list_paf_path,
    save_suffix,
    min_match_length, 
    min_idy
):
    with Parallel(n_jobs = len(list_paf_path)) as parallel:
        parallel(delayed(run_awk_for_filtering)(
            path_paf, path_paf.replace(".paf", save_suffix), min_match_length, min_idy
        )for path_paf in list_paf_path)

def run_awk_for_filtering(path_paf, path_save, min_match_length, min_idy):
    cmd = f"""awk '{{if(($10>{min_match_length}) && (($10/$11)>{min_idy})) print $0}}' {path_paf} > {path_save}"""
    subprocess.run(cmd, shell = True)

#%%
if __name__ == "__main__":
    hap = '1'
    min_match_length = 100
    min_idy = 0.5

    list_paf_name = [
        "Hifiasm_Binned_pri",
        "Hifiasm_Binned_alt",
        "Flye_HiFi_Binned",
        "Flye_ONT_Binned",
        "Unplaced_Scaffold",
        "Flye_HiFi_Binned_with0",
        "Flye_ONT_Binned_with0",
        "Flye_HiFi_Binned_with0a",
        "Flye_ONT_Binned_with0a",
        "Flye_ONT_Binned_short"
    ]
    list_paf_path = list(map(lambda name:
                             f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon/Round13.TGS_Gapcloser_Modified_V2_re.Unpl_Scf.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.Manual_Scaff.HiFi_trio.Hifiasm_Loc.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.ONT_trio.ONT_Amb.Manual_Scaff/Prealignment_for_Candidate_Filtering/Hap{hap}/Align.Option_Adjust/Align.REF_KOREFDraftASM_Hap{hap}_Split.QRY_{name}.asm20Pk15w5e0U1r50020knolongjoins100z800100qoccfrac0f00001.paf",
                             list_paf_name))
    # list_paf_path = list(map(lambda name:
    #                          f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon/Round13.TGS_Gapcloser_Modified_V2_re.Unpl_Scf.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.Manual_Scaff.HiFi_trio.Hifiasm_Loc.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.ONT_trio.ONT_Amb.Manual_Scaff/Prealignment_for_Candidate_Filtering/Hap{hap}/Align/Align.REF_KOREFDraftASM_Hap{hap}_Split.QRY_{name}.asm5P.paf",
    #                          list_paf_name))
    for path in list_paf_path:
        assert os.path.exists(path), path

    save_suffix = f".filt_match{min_match_length}_idy{min_idy}.paf"


    run(
        list_paf_path,
        save_suffix,
        min_match_length,
        min_idy
    )