#%%
import subprocess

PATH_BEDTOOLS = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/Bedtools/bedtools"
# %%
path_bedpe = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Repeat_Annotation/Segmental_Duplication/Biser/Hap1/FinalASM.KOREF_S1.Trio_Specific.hap1.Reordered.masked.SD.check_overlap.filt_overlap_50.err_lt30"
path_save = f"{path_bedpe}.merged.bed"
# %%
cmd = f"""awk 'BEGIN{{OFS="\\t"}} {{print $1, $2, $3; print $4, $5, $6}}' {path_bedpe} | {PATH_BEDTOOLS} sort -i - | {PATH_BEDTOOLS} merge -i - > {path_save}"""

subprocess.run(cmd, shell = True)
# %%
