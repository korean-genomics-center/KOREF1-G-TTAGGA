#%%
import gffutils

gff_file = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/Hap1.Blood_Cellline.CGP_PB.Add_Liftoff_CHM13_GRCh38/KOREF1_ASM_Hap1.CAT__CHM13_Liftoff__GRCh38_Liftoff.liftoff_final.gff3"

disable_transcripts = False
disable_genes = False

feature_db = gffutils.create_db(gff_file, gff_file + "_db", merge_strategy="create_unique", force=True,disable_infer_transcripts=disable_transcripts, disable_infer_genes=disable_genes, verbose=False)