#%%
import os,sys
from itertools import chain

from pathlib import Path
import pandas as pd

sys.path.append(str(Path(__file__).parents[1]))
from gff_parser import read_gff_as_dict

# Check novel transcripts are similar to annotated transcript/gene
# Check TPM value on cell-line or blood RNA-seq
# %%
def run(
        path_lowexpressed_transcripts,
        path_unknown_only_gff,
        dir_salmon_index,
        path_save_similar_to_known,
        path_save_lowexpressed_unknown
):
    os.makedirs(os.path.dirname(path_save_similar_to_known), exist_ok=True)
    os.makedirs(os.path.dirname(path_save_lowexpressed_unknown), exist_ok=True)

    list_lowexpressed_transcripts = read_textfile_as_list(path_lowexpressed_transcripts)
    dict_retained_to_duplicated_transcripts = get_deduplicated_status(dir_salmon_index)

    dict_unknown_gene2transcript = get_gene2transcript_from_gff(path_unknown_only_gff)
    list_unknown_transcripts = list(chain(*dict_unknown_gene2transcript.values()))

    list_unknown_transcripts_similar_to_known = get_transcripts_similar_to_anothers(dict_retained_to_duplicated_transcripts, list_unknown_transcripts)
    
    list_total_lowexpressed_transcripts = get_transcripts_extending_duplicates(dict_retained_to_duplicated_transcripts, list_lowexpressed_transcripts)
    list_unknown_lowexpressed_transcripts = list(set(list_total_lowexpressed_transcripts) & set(list_unknown_transcripts))

    save_list_values_into_file(list_unknown_transcripts_similar_to_known, path_save_similar_to_known)
    save_list_values_into_file(list_unknown_lowexpressed_transcripts, path_save_lowexpressed_unknown)


def read_textfile_as_list(path_file, strip = '\n'):
    list_values = list()
    with open(path_file, 'r') as fr:
        for line in fr:
            list_values.append(line.strip(strip))
    return list_values

def get_deduplicated_status(dir_salmon_index):
    path_salmon_duplicated_transcript = os.path.join(dir_salmon_index, "duplicate_clusters.tsv")

    table_duplicated_transcripts = pd.read_csv(path_salmon_duplicated_transcript, sep = '\t')
    dict_retained_to_duplicated = dict()
    for _, row in table_duplicated_transcripts.iterrows():
        retainedid = row["RetainedRef"]
        duplicateid = row["DuplicateRef"]
        if dict_retained_to_duplicated.get(retainedid) == None:
            dict_retained_to_duplicated[retainedid] = list()
        dict_retained_to_duplicated[retainedid].append(duplicateid)
    return dict_retained_to_duplicated

def get_gene2transcript_from_gff(path_gff):
    dict_gff = read_gff_as_dict(path_gff)

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
    return dict_gene_to_transcript

def get_transcripts_similar_to_anothers(dict_duplicates, list_transcripts):
    list_transcripts_similar_to_other = list()

    for retained, list_dedup in dict_duplicates.items():
        list_dupset = list_dedup + [retained]

        dupset_in_given_list = list(set(list_dupset) & set(list_transcripts))
        if len(dupset_in_given_list) > 0:
            if len(set(list_dupset) - set(list_transcripts)) > 0:
                list_transcripts_similar_to_other.extend(dupset_in_given_list)
    return list_transcripts_similar_to_other

def get_transcripts_extending_duplicates(dict_duplicates, list_transcripts):
    list_transcripts_extended = list()

    for trans_id in list_transcripts:
        list_transcripts_extended.append(trans_id)
        if dict_duplicates.get(trans_id) != None:
            list_transcripts_extended.extend(dict_duplicates[trans_id])
    return list_transcripts_extended

def save_list_values_into_file(list_values, path_save):
    with open(path_save, 'w') as fw:
        for val in list_values:
            fw.write(val + '\n')


if __name__ == "__main__":
    tpm_cutoff = 1

    for hap in [1, 2]:
        path_lowexpressed_transcripts = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/CDS_Fasta/Hap{hap}.Blood_Cellline.CGP_PB/Salmon/Quant_Salmon/Filter_by_TPM/Transcript_Filtering.Hap{hap}.TPM_cutoff_{tpm_cutoff}.transcripts_notexpressed.txt"

        path_unknown_only_gff = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/Process_Unknown/Hap{hap}.Blood_Cellline.CGP_PB/KOREF1_ASM_Hap{hap}.unknown_full.gff3"

        dir_salmon_index = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/CDS_Fasta/Hap{hap}.Blood_Cellline.CGP_PB/Salmon/Index_Salmon/KOREF1_ASM_Hap{hap}.whole_gff.salmon_index"

        path_save_similar_to_known = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/CDS_Fasta/Hap{hap}.Blood_Cellline.CGP_PB/Salmon/Quant_Salmon/Filter_Unknown/Transcript_Filtering.unknown_transcripts.Hap{hap}.similar_to_known.txt"
        path_save_lowexpressed_unknown = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/CDS_Fasta/Hap{hap}.Blood_Cellline.CGP_PB/Salmon/Quant_Salmon/Filter_Unknown/Transcript_Filtering.unknown_transcripts.Hap{hap}.TPM_cutoff_{tpm_cutoff}.lowexpression.txt"

        run(
            path_lowexpressed_transcripts,
            path_unknown_only_gff,
            dir_salmon_index,
            path_save_similar_to_known,
            path_save_lowexpressed_unknown
        )