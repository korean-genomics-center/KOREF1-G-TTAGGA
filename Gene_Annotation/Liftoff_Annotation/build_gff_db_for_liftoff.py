#%%
import gffutils

gff_file = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/CHM13/Gene_Annotation/chm13v2.0_RefSeq_Liftoff_v5.2.agat.gtf"

disable_transcripts = False
disable_genes = False

feature_db = gffutils.create_db(gff_file, gff_file + "_db", merge_strategy="create_unique", force=True,disable_infer_transcripts=disable_transcripts, disable_infer_genes=disable_genes, verbose=True)