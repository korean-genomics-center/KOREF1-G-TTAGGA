#%%
import gffutils

gff_file = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/CHM13/Gene_Annotation/chm13.draft_v2.0.gene_annotation.gff3"

disable_transcripts = False
disable_genes = False

feature_db = gffutils.create_db(gff_file, gff_file + "_db", merge_strategy="create_unique", force=True,disable_infer_transcripts=disable_transcripts, disable_infer_genes=disable_genes, verbose=False)