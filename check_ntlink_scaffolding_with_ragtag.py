#%%
import pandas as pd


def read_agp_file(path_file):
    def check_metainfo_lines(path_file):
        cnt_meta = 0
        with open(path_file, 'r') as fr:
            for line in fr:
                if line.startswith('#'):
                    cnt_meta += 1
                else:
                    break
        return cnt_meta
    
    lines_meta = check_metainfo_lines(path_file)

    table_agp = pd.read_csv(path_file, skiprows = lines_meta, sep = '\t', names = ["obj", "obj_beg", "obj_end", "part_number", "comp_type", "comp_id", "comp_beg", "comp_end", "ori"])
    return table_agp

def get_linked_components_from_agp(path_agp):
    table_agp = read_agp_file(path_agp)
    table_agp_nongap = table_agp[~table_agp["comp_type"].isin(['U', 'N'])]

    dict_linkage = dict()

    for obj in table_agp_nongap["obj"].unique():
        table_agp_nongap_of_obj = table_agp_nongap[table_agp_nongap["obj"] == obj]
        list_linked_comps = table_agp_nongap_of_obj["comp_id"].to_list()
        dict_linkage[obj] = list_linked_comps
    return dict_linkage

def relink_linkage_after_ntlink_round(dict_linkage_bef, dict_linkage_new):
    dict_linkage = dict()
    for obj, list_comp in dict_linkage_new.items():
        list_comp_update = list()
        for comp in list_comp:
            list_comp_update.extend(dict_linkage_bef.get(comp, [comp]))
        dict_linkage[obj] = list_comp_update
    return dict_linkage

def get_ntlink_linkage(path_ntlink_agp_format, n_rounds = 5):
    dict_linkage = dict()

    for r in range(1, n_rounds+1):
        path_ntlink_agp_round = path_ntlink_agp_format.format(round = '.'.join(["ntLink"]*r))

        dict_linkage_round = get_linked_components_from_agp(path_ntlink_agp_round)
        dict_linkage = relink_linkage_after_ntlink_round(dict_linkage, dict_linkage_round)
    return dict_linkage

#%%
path_ntlink_agp_format = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/ntLink.ONT_triobin_specific.round5_gapfilling/KOREF_S1.ul_L100kb_QV15.hifi_L1000_QV20.triobinning.asm.dip.hap1.p_ctg.fa.k32.w100.z1000.{round}.gap_fill.fa.agp"
n_rounds = 5

path_ragtag_agp = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/RagTag.default_option.CHM13/DraftASM.KOREF_S1.ul_L100kb_QV15.default_option.hap1.RagTag.default_option.CHM13/ragtag.scaffold.agp"
path_ragtag_agp_with_ntlink = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Scaffolding/Hifiasm/KOREF_S1.ul_L100kb_QV15.default_option/ntLink.ONT_triobin_specific.round5_gapfilling.RagTag.default_option.CHM13/DraftASM.KOREF_S1.ul_L100kb_QV15.default_option.hap1.ntLink.ONT_triobin_specific.round5_gapfilling.RagTag.default_option.CHM13/ragtag.scaffold.agp"
# %%
dict_final_ntlink_linkage = get_ntlink_linkage(path_ntlink_agp_format, n_rounds)
dict_ragtag_linkage_with_ntlink = get_linked_components_from_agp(path_ragtag_agp_with_ntlink)
dict_ragtag_linkage = get_linked_components_from_agp(path_ragtag_agp)
# %%
list_contigs_not_chr = list(filter(lambda name: not name.startswith("chr"), dict_ragtag_linkage.keys()))
# %%
for linkname, list_comp in dict_final_ntlink_linkage.items():
    if not linkname.startswith("ntLink"):
        continue
    for contig in list_contigs_not_chr:
        if contig in list_comp:
            print(linkname, contig)
# %%
