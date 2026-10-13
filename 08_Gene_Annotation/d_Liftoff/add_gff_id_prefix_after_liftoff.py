#%%
import os,sys
from copy import deepcopy

from pathlib import Path

sys.path.append(os.path.join(str(Path(__file__).parents[1]), "CAT_Pipeline", "Postprocess_CAT"))
from gff_parser import read_gff_as_dict, save_dict_gff_into_file

#%%
def run(path_gff, prefix, path_save_gff):
    os.makedirs(os.path.dirname(path_save_gff), exist_ok = True)
    dict_gff = read_gff_as_dict(path_gff)
    dict_gff_mod = add_prefix_to_dict_gff(dict_gff, prefix)
    save_dict_gff_into_file(dict_gff_mod, path_save_gff)

def add_prefix_to_dict_gff(dict_gff, prefix):
    dict_gff_mod = {
        "metainfo":deepcopy(dict_gff["metainfo"]),
        "data":list()
    }

    for dict_data in dict_gff["data"]:
        dict_data_mod = deepcopy(dict_data)
        dict_data_mod["attrib"]["ID"] = f"{prefix}:{dict_data['attrib']['ID']}"
        if dict_data["attrib"].get("Parent") != None:
            dict_data_mod["attrib"]["Parent"] = f"{prefix}:{dict_data['attrib']['Parent']}"
        dict_gff_mod["data"].append(dict_data_mod)
    return dict_gff_mod


#%%

if __name__ == "__main2__":
    path_gff = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/Liftoff/From.GRCh38_primary_GENCODEv49.To.KOREF1_G_ASM_Hap2/liftoff_annotation.KOREF1_G_ASM_Hap2.from_GRCh38_primary_GENCODEv49.annotate.gff3_polished"
    prefix = "GRCh38"
    path_save = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/Liftoff/From.GRCh38_primary_GENCODEv49.To.KOREF1_G_ASM_Hap2/liftoff_annotation.KOREF1_G_ASM_Hap2.from_GRCh38_primary_GENCODEv49.annotate.prefix_add.gff3_polished"
    run(path_gff, prefix, path_save)
# %%

if __name__ == "__main__":
    for hap in [1, 2]:
        path_gff = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/Liftoff/From.CHM13v2_GENCODEv35_CATLiftoff.To.KOREF1_G_ASM_Hap{hap}/liftoff_annotation.KOREF1_G_ASM_Hap{hap}.from_CHM13v2_GENCODEv35_CATLiftoff.annotate.gff3_polished"
        prefix = "CHM13"
        path_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/Liftoff/From.CHM13v2_GENCODEv35_CATLiftoff.To.KOREF1_G_ASM_Hap{hap}/liftoff_annotation.KOREF1_G_ASM_Hap{hap}.from_CHM13v2_GENCODEv35_CATLiftoff.annotate.prefix_add.gff3_polished"
        run(path_gff, prefix, path_save)