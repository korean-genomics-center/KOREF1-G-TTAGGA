#%%
import pandas as pd

path_gff_chm13 = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/ReferenceGenome/CHM13/Gene_Annotation/chm13v2.0_RefSeq_Liftoff_v5.2.agat.gtf"
path_gff_hap1 = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/Liftoff/From.chm13v2_RefSeq_Liftoff_v5.2.To.KOREF_ASM_20251027_Hap1/liftoff_annotation.KOREF_ASM_Hap1.from_chm13v2_refseq_v5.2.annotate.gtf_polished"
path_gff_hap2 = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/Liftoff/From.chm13v2_RefSeq_Liftoff_v5.2.To.KOREF_ASM_20251027_Hap2/liftoff_annotation.KOREF_ASM_Hap2.from_chm13v2_refseq_v5.2.annotate.gtf_polished"
# %%
table_gff_chm13 = pd.read_csv(path_gff_chm13, sep = '\t', names = ["chr", "source", "feature", "start", "end", "score", "strand", "frame", "attrib"], skiprows = 2)
table_gff_hap1 = pd.read_csv(path_gff_hap1, sep = '\t', names = ["chr", "source", "feature", "start", "end", "score", "strand", "frame", "attrib"], skiprows = 2)
table_gff_hap2 = pd.read_csv(path_gff_hap2, sep = '\t', names = ["chr", "source", "feature", "start", "end", "score", "strand", "frame", "attrib"], skiprows = 2)
# %%
def get_attrib_key_from_value(value, key):
    list_values = value.strip().strip(';').split('";')
    for val in list_values:
        val = val.strip()
        val_key = val.split()[0]

        if val_key == key:
            val_key, val_val = val.split(' "')
            val_key = val_key.strip()
            val_val = val_val.strip('"')
            return val_val
    return None
#%%
table_gff_chm13_gene = table_gff_chm13[table_gff_chm13["feature"] == "gene"]
table_gff_hap1_gene = table_gff_hap1[table_gff_hap1["feature"] == "gene"]
table_gff_hap2_gene = table_gff_hap2[table_gff_hap2["feature"] == "gene"]
# %%
table_gff_chm13_gene["biotype"] = table_gff_chm13_gene["attrib"].apply(lambda val: get_attrib_key_from_value(val, "gene_biotype"))
# %%

table_gff_hap1_gene["biotype"] = table_gff_hap1_gene["attrib"].apply(lambda val: get_attrib_key_from_value(val, "gene_biotype"))
table_gff_hap2_gene["biotype"] = table_gff_hap2_gene["attrib"].apply(lambda val: get_attrib_key_from_value(val, "gene_biotype"))
# %%
