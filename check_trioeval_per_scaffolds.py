#%%
import pandas as pd

#%%
# path_trioeval_s = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/TrioEval/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.asm.dip.hap1.p_ctg.parental_specific_kmer.trioeval_err_stdout.S_line"
path_trioeval_s = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option/TrioEval/KKOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.Trio_Specific.asm.dip.hap1.p_ctg.parental_only_kmer.trioeval_err_stdout.S_line"

col_kmercount_hap = "#patKmer"
table_trioeval = pd.read_csv(path_trioeval_s, delim_whitespace = True, names = ["S",  "seqName",     "#patKmer",  "#matKmer",  "#pat-pat",  "#pat-mat",  "#mat-pat",  "#mat-mat",  "seqLen"])
# %%
table_trioeval["totKmer"] = table_trioeval["#patKmer"] + table_trioeval["#matKmer"]
table_trioeval["totPair"] = table_trioeval["#pat-pat"] + table_trioeval["#pat-mat"] + table_trioeval["#mat-pat"] + table_trioeval["#mat-mat"]
# %%
table_trioeval["errRate"] = table_trioeval.apply(lambda row: None if row["totKmer"] == 0 else 1 - row[col_kmercount_hap] / row["totKmer"], axis = 1)
table_trioeval["switchErrRate"] = table_trioeval.apply(lambda row: None if row["totPair"] == 0 else (row["#pat-mat"] + row["#mat-pat"]) / row["totPair"], axis = 1)
# %%
from matplotlib import pyplot as plt

table_trioeval_nonnaerr = table_trioeval.dropna(subset = ["errRate"])
table_trioeval_nonnaswerr = table_trioeval.dropna(subset = ["switchErrRate"])

plt.scatter(
    table_trioeval_nonnaerr["seqLen"], table_trioeval_nonnaerr["errRate"]
)
# %%
plt.scatter(
    table_trioeval_nonnaswerr["seqLen"], table_trioeval_nonnaswerr["switchErrRate"]
)
# %%
