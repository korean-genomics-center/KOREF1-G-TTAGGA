#%%
import os
import pandas as pd
import numpy as np

# 1. Keep all expression-quantified annotations
# 2. Remove spurious protein predictions
# 3. If there’s no Diamond/InterProScan hit, remove too short transcripts
# 4. If CDS length is over 100bp but there’s no Diamond/InterProScan hit, report it as a putative false-positive 

def run(path_evaluation_summary, path_output_prefix):
    os.makedirs(os.path.dirname(path_output_prefix), exist_ok = True)

    table_evaluation_summary = pd.read_csv(path_evaluation_summary, sep = '\t')

    table_noexpquant = table_evaluation_summary[~table_evaluation_summary["Expression_Quant"]]

    table_spurious = table_noexpquant[~table_noexpquant["NonSpurious"]]

    table_nohitshort = table_noexpquant[
        np.logical_and(
            np.logical_and(
                ~table_noexpquant["Diamond"],
                ~table_noexpquant["InterProScan"]
            ),
            ~table_noexpquant["CDS_Length"]
        )
    ]
    table_nohitlong = table_noexpquant[
        np.logical_and(
            np.logical_and(
                np.logical_and(
                    ~table_noexpquant["Diamond"],
                    ~table_noexpquant["InterProScan"]
                ),
                table_noexpquant["CDS_Length"]
            ),
            table_noexpquant["NonSpurious"]
        )
    ]

    list_transcripts_filtout = list(set(table_spurious["Transcript_ID"]) | set(table_nohitshort["Transcript_ID"]))
    list_transcripts_warn = table_nohitlong["Transcript_ID"].to_list()

    print("Filtout", len(list_transcripts_filtout))
    print("Warn", len(list_transcripts_warn))

    write_list_values_as_txt(list_transcripts_filtout, f"{path_output_prefix}.filtout.txt")
    write_list_values_as_txt(list_transcripts_warn, f"{path_output_prefix}.warn.txt")

def write_list_values_as_txt(list_values, path_save):
    with open(path_save, 'w') as fw:
        for val in list_values:
            fw.write(f"{val}\n")

if __name__ == "__main__":
    for hap in [1, 2]:
        path_evaluation_summary = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/Process_Unknown/Hap{hap}.Blood_Cellline.CGP_PB/Evaluate.KOREF1_ASM_Hap{hap}.unknown_transcript.TPM_cutoff_1.tsv"
        path_output_prefix = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/Process_Unknown/Hap{hap}.Blood_Cellline.CGP_PB/Evaluate.KOREF1_ASM_Hap{hap}.unknown_transcript.TPM_cutoff_1.transcrpts"

        run(path_evaluation_summary, path_output_prefix)