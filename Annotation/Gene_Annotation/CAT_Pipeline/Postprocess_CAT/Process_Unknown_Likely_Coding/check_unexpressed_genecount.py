#%%
import os,sys

from pathlib import Path
import pandas as pd

sys.path.append(str(Path(__file__).parents[1]))
from gff_parser import read_gff_as_dict

#%%
dict_gff = read_gff_as_dict("/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Run_CAT/Hap1.Blood_Cellline.CGP_PB/consensus_gene_set/KOREF1_ASM_Hap1.gff3")

#%%
list_transcripts_lowexpression = list()
with open("/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/CDS_Fasta/Hap1.Blood_Cellline.CGP_PB/Salmon/Quant_Salmon/Filter_by_TPM/Transcript_Filtering.Hap1.TPM_cutoff_1.transcripts_notexpressed.txt", 'r') as fr:
    for line in fr:
        list_transcripts_lowexpression.append(line.strip('\n'))

#%%
table_duplicated_transcripts = pd.read_csv("/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/CDS_Fasta/Hap1.Blood_Cellline.CGP_PB/Salmon/Index_Salmon/KOREF1_ASM_Hap1.whole_gff.salmon_index/duplicate_clusters.tsv", sep = '\t')
dict_retained_to_duplicated = dict()
for _, row in table_duplicated_transcripts.iterrows():
    retainedid = row["RetainedRef"]
    duplicateid = row["DuplicateRef"]
    if dict_retained_to_duplicated.get(retainedid) == None:
        dict_retained_to_duplicated[retainedid] = list()
    dict_retained_to_duplicated[retainedid].append(duplicateid)

#%%
list_transcripts_lowexpression_extend = list()
for transcid in list_transcripts_lowexpression:
    list_transcripts_lowexpression_extend.append(transcid)
    if dict_retained_to_duplicated.get(transcid) != None:
        list_transcripts_lowexpression_extend.extend(dict_retained_to_duplicated[transcid])
#%%
list_protein_coding_gene = list()
list_unknown_gene = list()

for gff_line_info in dict_gff["data"]:
    if gff_line_info["feat_type"] == "gene":
        gene_id = gff_line_info["attrib"]["ID"]
        gene_type = gff_line_info["attrib"]["gene_biotype"]

        if gene_type == "protein_coding":
            list_protein_coding_gene.append(gene_id)
        if gene_type == "unknown_likely_coding":
            list_unknown_gene.append(gene_id)
#%%
dict_gene_to_transcript = dict()

for gff_line_info in dict_gff["data"]:
    if gff_line_info["feat_type"] == "gene":
        gene_id = gff_line_info["attrib"]["ID"]
        assert dict_gene_to_transcript.get(gene_id) == None, gene_id
        dict_gene_to_transcript[gene_id] = list()
    elif gff_line_info["feat_type"] == "transcript":
        parent_of_feat = gff_line_info["attrib"]["Parent"]
        transc_id = gff_line_info["attrib"]["ID"]

        assert dict_gene_to_transcript.get(parent_of_feat) != None, gene_id
        dict_gene_to_transcript[parent_of_feat].append(transc_id)

#%%
set_lowexpressed_transcripts = set(list_transcripts_lowexpression_extend)

list_lowexpressed_genes = list()
for geneid, list_transcrips in dict_gene_to_transcript.items():
    if set(list_transcrips) <= set_lowexpressed_transcripts:
        list_lowexpressed_genes.append(geneid)
# %%
print(len(list_lowexpressed_genes))
print(len(list_lowexpressed_genes)/len(dict_gene_to_transcript)*100)

# %%
print(len(set(list_lowexpressed_genes)&set(list_protein_coding_gene)))
print(len(set(list_lowexpressed_genes)&set(list_protein_coding_gene))/len(list_protein_coding_gene)*100)
# %%
print(len(set(list_lowexpressed_genes)&set(list_unknown_gene)))
print(len(set(list_lowexpressed_genes)&set(list_unknown_gene))/len(list_unknown_gene)*100)

# %%
