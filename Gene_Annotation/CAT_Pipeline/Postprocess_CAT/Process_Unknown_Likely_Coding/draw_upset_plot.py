#%%
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import upsetplot
from upsetplot import from_indicators
from upsetplot import UpSet
#%%
tpm_cutoff = 1
path_evaluation_result = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/Process_Unknown/Hap1.Blood_Cellline.CGP_PB/Evaluate.KOREF1_ASM_Hap1.unknown_transcript.TPM_cutoff_{tpm_cutoff}.tsv"

table_evaluation_result = pd.read_csv(path_evaluation_result, sep = '\t')

#%%
list_indicators = list(table_evaluation_result.columns)[1:]

plt.rcParams["font.size"] = 13
UpSet(
    from_indicators(
        list_indicators, data=table_evaluation_result
    ), show_counts=True
)
# %%
