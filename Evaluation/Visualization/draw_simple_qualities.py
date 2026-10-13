#%%
from matplotlib import pyplot as plt


#%%
dict_qv = {
    "KOREF S1 (Paternal)": 81.1933,
    "KOREF S1 (Maternal)": 79.0262,
    "CHM13 v1.1": 73.94,
    "CN1 (Paternal)": 59.36,
    "CN1 (Maternal)": 60.10,
    "T2T-YAO (Paternal)": 72.28,
    "T2T-YAO (Maternal)": 70.49,
    "I002C (Paternal)": 83.08,
    "I002C (Maternal)": 82.05
}


plt.rcParams["font.size"] = 15
plt.figure(figsize = (7, 3))
plt.bar(
    dict_qv.keys(),
    dict_qv.values(),
    color = "gray"
)
plt.ylabel("QV")
plt.xticks(rotation = -75, ha = "center", fontsize = plt.rcParams["font.size"]-1)

for ind, val in enumerate(dict_qv.values()):
    plt.annotate(f"{val:.2f}", (ind, val), va = "bottom", ha = "center", fontsize = plt.rcParams["font.size"]-2)
    
plt.ylim(top = 90)


plt.show()
# %%

dict_gci = {
    "KOREF1-G-TTAGGA v1 (Paternal)": 84.9181,
    "KOREF1-G-TTAGGA v1 (Maternal)": 80.7832,
    "CHM13 v2": 87.04,
    "CN1 v1.0.1 (Paternal)": 77.90,
    "CN1 v1.0.1 (Maternal)": 66.79
}

plt.rcParams["font.size"] = 15
plt.figure(figsize = (5, 3))
plt.barh(
    list(dict_gci.keys())[::-1],
    list(dict_gci.values())[::-1],
    color = "gray"
)
plt.xlabel("GCI score")
# plt.yticks(rotation = -75, ha = "center", fontsize = plt.rcParams["font.size"]-1)

for ind, val in enumerate(list(dict_gci.values())[::-1]):
    plt.annotate(f"{val:.2f}", (val+1, ind), va = "center", ha = "left", fontsize = plt.rcParams["font.size"]-1)
plt.xlim(right = 105)
# plt.ylim(top = 95)

plt.savefig("fig_gci.png", dpi = 350, bbox_inches = "tight")
plt.show()
# %%
