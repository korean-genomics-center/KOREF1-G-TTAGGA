#%%
from matplotlib_venn import venn3_unweighted
from matplotlib import pyplot as plt

def read_paf(path_paf):
    table_paf = pd.read_csv(path_paf, sep = '\t', names = list(range(1, 13)), usecols = list(range(1, 13)), low_memory = False)
    table_paf = table_paf.rename(
        columns = {
            1:"qryname",
            2:"qrylen",
            3:"qry_alignstart",
            4:"qry_alignend",
            5:"strand",
            6:"refname",
            7:"reflen",
            8:"ref_alignstart",
            9:"ref_alignend",
            10:"match",
            11:"align",
            12:"quality"
        }
    )
    return table_paf

def draw_venn_with_filter(table1, table2, table3, name1, name2, name3, ref, idy = 0, front = True, ratio = 0.3):
    table1_ref = table1[table1["refname"] == ref]
    table2_ref = table2[table2["refname"] == ref]
    table3_ref = table3[table3["refname"] == ref]

    reflen = table1_ref["reflen"].to_list()[0]
    if front:
        table1_ref = table1_ref[table1_ref["ref_alignstart"] < reflen*ratio]
        table2_ref = table2_ref[table2_ref["ref_alignstart"] < reflen*ratio]
        table3_ref = table3_ref[table3_ref["ref_alignstart"] < reflen*ratio]
    else:
        table1_ref = table1_ref[(table1_ref["reflen"]-table1_ref["ref_alignend"]) < reflen*ratio]
        table2_ref = table2_ref[(table2_ref["reflen"]-table2_ref["ref_alignend"]) < reflen*ratio]
        table3_ref = table3_ref[(table3_ref["reflen"]-table3_ref["ref_alignend"]) < reflen*ratio]

    venn3_unweighted(
        [
            set(table1_ref[table1_ref["identity"] > idy]["qryname"]),
            set(table2_ref[table2_ref["identity"] > idy]["qryname"]),
            set(table3_ref[table3_ref["identity"] > idy]["qryname"])
        ],
        set_labels = (name1, name2, name3)
    )
    plt.show()



#%%
table_paf_ori = pd.read_csv("/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon/Round13.TGS_Gapcloser_Modified_V2_re.Unpl_Scf.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.Manual_Scaff.HiFi_trio.Hifiasm_Loc.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.ONT_trio.ONT_Amb.Manual_Scaff/Prealignment_for_Candidate_Filtering/Hap1/Align/Align.REF_KOREFDraftASM_Hap1_Split.QRY_Any.asm5P.filtered.idy_0.7.align_5000.align2hang_0.5.paf", sep = '\t')
table_paf_noasm = pd.read_csv("/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon/Round13.TGS_Gapcloser_Modified_V2_re.Unpl_Scf.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.Manual_Scaff.HiFi_trio.Hifiasm_Loc.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.ONT_trio.ONT_Amb.Manual_Scaff/Prealignment_for_Candidate_Filtering/Hap1/Align.Option_Adjust/Align.REF_KOREFDraftASM_Hap1_Split.QRY_Any.asm5P.filtered.idy_0.7.align_5000.align2hang_0.5.paf", sep = '\t')
table_paf_asm20 = pd.read_csv("/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gapfilling/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TGS_Gapcloser_wo_Racon/Round13.TGS_Gapcloser_Modified_V2_re.Unpl_Scf.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.Manual_Scaff.HiFi_trio.Hifiasm_Loc.Hifiasm_Bin.Flye_HiFi_Bin.Flye_ONT_Bin.ONT_trio.ONT_Amb.Manual_Scaff/Prealignment_for_Candidate_Filtering/Hap1/Align.Option_Adjust/Align.REF_KOREFDraftASM_Hap1_Split.QRY_Any.asm20Pk15w5e0U1r50020knolongjoin.filtered.idy_0.7.align_5000.align2hang_0.5.paf", sep = '\t')

#%%
draw_venn_with_filter(table_paf_ori, table_paf_noasm, table_paf_asm20, "ori", "noasm", "asm20", 26, idy = 0.85, front = False)
draw_venn_with_filter(table_paf_ori, table_paf_noasm, table_paf_asm20, "ori", "noasm", "asm20", 28, idy = 0.85, front = True)