#%%
import os,sys
from itertools import chain

from pathlib import Path
import pandas as pd

sys.path.append(str(Path(__file__).parents[1]))
from gff_parser import read_gff_as_dict


#%%
def run(
    path_unknown_transcript_gff,
    path_diamond_result,
    path_interproscan_result,
    path_interproscan_result_antifam,
    path_isoseq_supp_gff,
    path_check_duplicated_with_known,
    path_lowexpression_check,
    path_shortcds_check,
    path_save_total_check
):
    os.makedirs(os.path.dirname(path_save_total_check), exist_ok = True)
    list_total_unknown_transcripts = get_transcript_ids_from_gff(path_unknown_transcript_gff)

    list_diamond_searched_transcripts = get_transcript_ids_from_table_with_col_ind(path_diamond_result, 0)
    list_interproscan_searched_transcripts = get_transcript_ids_from_table_with_col_ind(path_interproscan_result, 0)
    list_spurious_transcripts = get_transcript_ids_from_table_with_col_ind(path_interproscan_result_antifam, 0)
    list_isoseq_supported_transcripts = get_transcript_ids_from_gff(path_isoseq_supp_gff)
    list_duplicated_known_transcripts = get_transcript_ids_from_table_with_col_ind(path_check_duplicated_with_known, 0)
    list_lowexpressed_transcripts = get_transcript_ids_from_table_with_col_ind(path_lowexpression_check, 0)
    list_tooshort_transcripts = get_transcript_ids_from_table_with_col_ind(path_shortcds_check, 0)

    dict_evaluation_results = {
        "Diamond" : {"list":list_diamond_searched_transcripts, "include":True},
        "InterProScan" : {"list":list_interproscan_searched_transcripts, "include":True},
        "NonSpurious":{"list":list_spurious_transcripts, "include":False},
        "Isoseq_Support" : {"list":list_isoseq_supported_transcripts, "include":True},
        "Homology_Known" : {"list":list_duplicated_known_transcripts, "include":True},
        "Expression_Quant" : {"list":list_lowexpressed_transcripts, "include":False},
        "CDS_Length" : {"list":list_tooshort_transcripts, "include":False},
    }
    table_evaluation_result = make_table_with_filtering_evaluations(list_total_unknown_transcripts, dict_evaluation_results)
    table_evaluation_result.to_csv(path_save_total_check, sep = '\t', index = False)


def get_transcript_ids_from_table_with_col_ind(path_table, ind_col, sep = '\t'):
    with open(path_table, 'r') as fr:
        firstline = fr.readline().strip('\n')
        n_cols = len(firstline.split(sep))
    
    table = pd.read_csv(path_table, sep = sep, names = list(range(n_cols)))
    list_transcript_ids = table.loc[:, ind_col].to_list()
    return list_transcript_ids

def get_transcript_ids_from_gff(path_gff):
    dict_gff = read_gff_as_dict(path_gff)

    list_transcript_ids = list()
    for dict_data in dict_gff["data"]:
        if dict_data["feat_type"] == "transcript":
            list_transcript_ids.append(dict_data["attrib"]["ID"])
    return list_transcript_ids

def make_table_with_filtering_evaluations(list_total_transcripts, dict_evaluation_results):
    table = pd.DataFrame({"Transcript_ID":list_total_transcripts})

    for eval_name in dict_evaluation_results.keys():
        set_eval_transcripts = set(dict_evaluation_results[eval_name]["list"])

        if dict_evaluation_results[eval_name]["include"]:
            table[eval_name] = table["Transcript_ID"].apply(lambda val: val in set_eval_transcripts)
        else:
            table[eval_name] = table["Transcript_ID"].apply(lambda val: val not in set_eval_transcripts)
    return table


if __name__ == "__main__":
    tpm_cutoff = 0.5
    for hap in [1, 2]:
        path_unknown_transcript_gff = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/Process_Unknown/Hap{hap}.Blood_Cellline.CGP_PB/KOREF1_ASM_Hap{hap}.unknown_transcript.gff3"

        path_diamond_result = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/CDS_Fasta/Hap{hap}.Blood_Cellline.CGP_PB/Diamond_Blastp/Diamond_Blastp.KOREF1_ASM_Hap{hap}.unknown_gff.protein.against_swissprot.evalue_1e5.tsv"
        path_interproscan_result = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/CDS_Fasta/Hap{hap}.Blood_Cellline.CGP_PB/InterProScan/InterProScan_Result.KOREF1_ASM_Hap{hap}.unknown_gff.protein.db_without_Coils_MobiDBLite.Without_AntiFam.tsv"
        path_interproscan_result_antifam = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/CDS_Fasta/Hap{hap}.Blood_Cellline.CGP_PB/InterProScan/InterProScan_Result.KOREF1_ASM_Hap{hap}.unknown_gff.protein.db_without_Coils_MobiDBLite.Spurious_AntiFam.tsv"
        path_isoseq_supp_gff = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/Process_Unknown/Hap{hap}.Blood_Cellline.CGP_PB/KOREF1_ASM_Hap{hap}.unknown_transcript.isoseq_supported.gff3"
        path_check_duplicated_with_known = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/CDS_Fasta/Hap{hap}.Blood_Cellline.CGP_PB/Salmon/Quant_Salmon/Filter_Unknown/Transcript_Filtering.unknown_transcripts.Hap{hap}.similar_to_known.txt"
        path_lowexpression_check = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/CDS_Fasta/Hap{hap}.Blood_Cellline.CGP_PB/Salmon/Quant_Salmon/Filter_Unknown/Transcript_Filtering.unknown_transcripts.Hap{hap}.TPM_cutoff_{tpm_cutoff}.lowexpression.txt"
        path_shortcds_check = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/CDS_Fasta/Hap{hap}.Blood_Cellline.CGP_PB/CDS_Length/GFF_unknown/Transcript_Filtering.Hap{hap}.GFF_unknown.CDSLength_cutoff_100.tooshort.txt"

        path_save_total_check = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/Process_Unknown/Hap{hap}.Blood_Cellline.CGP_PB/Evaluate.KOREF1_ASM_Hap{hap}.unknown_transcript.TPM_cutoff_{tpm_cutoff}.tsv"

        run(
            path_unknown_transcript_gff,
            path_diamond_result,
            path_interproscan_result,
            path_interproscan_result_antifam,
            path_isoseq_supp_gff,
            path_check_duplicated_with_known,
            path_lowexpression_check,
            path_shortcds_check,
            path_save_total_check
        )