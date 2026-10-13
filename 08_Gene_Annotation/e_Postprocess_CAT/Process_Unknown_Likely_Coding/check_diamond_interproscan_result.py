#%%
import pandas as pd

path_diamond_result = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/CDS_Fasta/Hap2.Blood_Cellline.CGP_PB/Diamond_Blastp/Diamond_Blastp.KOREF1_ASM_Hap2.unknown_gff.protein.against_swissprot.evalue_1e5.tsv"
path_interproscan_result = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/CDS_Fasta/Hap2.Blood_Cellline.CGP_PB/InterProScan/InterProScan_Result.KOREF1_ASM_Hap2.unknown_gff.protein.db_without_Coils_MobiDBLite.Without_AntiFam.tsv"
path_interproscan_result_antifam = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/CDS_Fasta/Hap2.Blood_Cellline.CGP_PB/InterProScan/InterProScan_Result.KOREF1_ASM_Hap2.unknown_gff.protein.db_without_Coils_MobiDBLite.Spurious_AntiFam.tsv"

list_diamond_columns = ["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","qcovhsp","scovhsp","stitle"]
list_interproscan_columns = ["name", "seq_md5", "seq_len", "analysis", "signature_accession", "signature_description", "start", "end", "score", "status", "date", "interpro_accession", "interpro_description", "go", "pathway"]
# %%
table_diamond = pd.read_csv(path_diamond_result, sep = '\t', names = list_diamond_columns)
table_interproscan = pd.read_csv(path_interproscan_result, sep = '\t', names = list_interproscan_columns)
table_interproscan_antifam = pd.read_csv(path_interproscan_result_antifam, sep = '\t', names = list_interproscan_columns)
# %%
list_transcripts_passed_diamond = list(set(table_diamond["qseqid"]))
list_transcripts_passed_interproscan = list(set(table_interproscan["name"]))
list_transcripts_spurious = list(set(table_interproscan_antifam["name"]))

# %%
len(list_transcripts_passed_diamond)
# %%
len(list_transcripts_passed_interproscan)

# %%
len(set(list_transcripts_passed_diamond) & set(list_transcripts_passed_interproscan))
# %%
len(set(list_transcripts_passed_diamond) | set(list_transcripts_passed_interproscan))
# %%
from matplotlib_venn import venn3_unweighted
from matplotlib import pyplot as plt

fig = plt.figure(figsize = (7, 5))

venn3_unweighted(
    [set(list_transcripts_passed_diamond),
     set(list_transcripts_passed_interproscan),
     set(list_transcripts_spurious)],
     set_labels = ["Diamond", "InterProScan", "Spurious"]
)
plt.show()
# %%
