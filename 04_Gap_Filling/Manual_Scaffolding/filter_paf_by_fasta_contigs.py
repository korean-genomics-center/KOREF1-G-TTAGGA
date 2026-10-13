#%%
import os, sys
from joblib import Parallel, delayed
from pathlib import Path

sys.path.append(str(Path(__file__).parents[2] / "common"))

from fasta_parser import read_fasta_as_dict


def main(path_paf, path_fasta_ref, path_fasta_qry, path_save):
    dict_fasta_ref = read_fasta_as_dict(path_fasta_ref)
    dict_fasta_qry = read_fasta_as_dict(path_fasta_qry)

    set_seqname_ref = set(dict_fasta_ref.keys())
    set_seqname_qry = set(dict_fasta_qry.keys())

    with open(path_paf, 'r') as fr, open(path_save, 'w') as fw:
        for line in fr:
            list_values = line.split('\t')
            qryname = list_values[0]
            refname = list_values[5]

            if qryname in set_seqname_qry and refname in set_seqname_ref:
                fw.write(line)


if __name__ == "__main__":
    dict_hap_to_parental = {
        1:'p',
        2:'m'
    }

    qv_thresh = 40

    dict_binnedasm_format = {
        "Hifiasm_Binned_pri":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.Binned_Input/ONT_UL.Filter_L50kb_QV20_Q30Perc80.HiFi.Filter_L1000_QV20/Hap{hap}/KOREF_S1.Binned_Assembly.ul_L50kb_QV20_Q30Perc80.triobin_{triobin}.hifi_L1000_QV20.triobin_{triobin}.haploid.l0.asm.bp.p_ctg.qv{qv_thresh}.fasta",
        "Hifiasm_Binned_alt":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.Binned_Input/ONT_UL.Filter_L50kb_QV20_Q30Perc80.HiFi.Filter_L1000_QV20/Hap{hap}/KOREF_S1.Binned_Assembly.ul_L50kb_QV20_Q30Perc80.triobin_{triobin}.hifi_L1000_QV20.triobin_{triobin}.haploid.l0.asm.bp.a_ctg.qv{qv_thresh}.fasta",

        # "Hifiasm_Binned_woONT_pri":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20/Hap{hap}/KOREF_S1.Binned_Assembly.hifi_L1000_QV20.triobin_{triobin}.haploid.l0.asm.bp.p_ctg.qv{qv_thresh}.fasta",
        # "Hifiasm_Binned_woONT_alt":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20/Hap{hap}/KOREF_S1.Binned_Assembly.hifi_L1000_QV20.triobin_{triobin}.haploid.l0.asm.bp.a_ctg.qv{qv_thresh}.fasta",

        "Flye_HiFi_Binned":"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/HiFi.Filter_L1000_QV20/Hap{hap}/assembly.qv{qv_thresh}.fasta",
        "Flye_ONT_Binned" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Flye/KOREF_S1.Binned_Input/ONT_UL.Filter_L50kb_QV20_Q30Perc80/Hap{hap}/assembly.qv{qv_thresh}.fasta",

        "Unplaced_Scaffold" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/Chromosome_level_scaffolds.Assign_by_CHM13/CuratedDraftASM.KOREF_S1.ul_L100kb_QV15.Trio_Specific.ntLink.adjusted_telomere.hap{hap}.RagTag_with_Telomere_Assign.minimap2_latest.CHM13.unplaced.fasta"
    }

    path_paf_format = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon/Round13.TGS_Gapcloser_Modified_V2_re.Unpl_Scf.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.Manual_Scaff.HiFi_trio.Hifiasm_Loc.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.ONT_trio.ONT_Amb.Manual_Scaff/Prealignment_for_Candidate_Filtering/Hap{hap}/Align.Pairwise/Align.Hap{hap}.REF.{name_ref}.QRY.{name_qry}.asm5Pc_triobin_{triobin}.filt_match100.paf"

    path_save_format = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon_wo_RawONT/Contig_Alignment/Hap{hap}/Align.Pairwise/Filter_QV{qv_thresh}/Align.Hap{hap}.REF.{name1}.QRY.{name2}.asm5Pc_triobin_{triobin}.filt_match100.paf"

    list_option_comb = list()
    for hap in [1, 2]:
        triobin = dict_hap_to_parental[hap]
        
        list_binnedasm = list(dict_binnedasm_format.keys())
        for name1 in list_binnedasm:
            for name2 in list_binnedasm:
                list_option_comb.append(
                    (hap, triobin, name1, name2)
                )
    
    with Parallel(n_jobs = 20) as parallel:
        parallel(delayed(main)(
            path_paf_format.format(hap = hap, name_ref = name1, name_qry = name2, triobin = triobin),
            dict_binnedasm_format[name1].replace("{hap}", str(hap)).replace("{triobin}", triobin).replace("{qv_thresh}", str(qv_thresh)),
            dict_binnedasm_format[name2].replace("{hap}", str(hap)).replace("{triobin}", triobin).replace("{qv_thresh}", str(qv_thresh)),
            path_save_format.format(hap = hap, name1 = name1, name2 = name2, triobin = triobin, qv_thresh = qv_thresh)
        ) for (hap, triobin, name1, name2) in list_option_comb)
                
                # path_fasta1 = dict_binnedasm_format[name1].replace("{hap}", str(hap)).replace("{triobin}", triobin).replace("{qv_thresh}", str(qv_thresh))
                # path_fasta2 = dict_binnedasm_format[name2].replace("{hap}", str(hap)).replace("{triobin}", triobin).replace("{qv_thresh}", str(qv_thresh))

                # path_paf_input = path_paf_format.format(hap = hap, name_ref = name1, name_qry = name2, triobin = triobin)
                # path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon_wo_RawONT/Contig_Alignment/Hap{hap}/Align.Pairwise/Filter_QV{qv_thresh}/Align.Hap{hap}.REF.{name1}.QRY.{name2}.asm5Pc_triobin_{triobin}.filt_match100.paf"

                # main(path_paf_input, path_fasta1, path_fasta2, path_save)