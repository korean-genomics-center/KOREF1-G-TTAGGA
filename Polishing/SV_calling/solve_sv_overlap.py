#%%
from collections import Counter

import pandas as pd


COL_SAMPLE = "0_SAMPLE__1_SAMPLE"

def run(path_vcf, path_save, nsupp_without_vaf_comparison = 10):
    table_vcf = read_vcf_as_table(path_vcf)
    
    table_vcf["DV"] = table_vcf[COL_SAMPLE].apply(lambda val: int(val.split(':')[3]))
    table_vcf["DR"] = table_vcf[COL_SAMPLE].apply(lambda val: int(val.split(':')[4]))
    table_vcf["SVType"] = table_vcf[COL_SAMPLE].apply(lambda val: val.split(':')[2])
    
    table_vcf["SVLEN"] = table_vcf.apply(lambda row: len(row["ALT"]) - len(row["REF"]), axis = 1)
    
    list_varids_filtout = list()
    
    for chrname in table_vcf["#CHROM"].unique():
        list_filtout_on_chr = check_variant_overlap_on_single_chromosome(table_vcf[table_vcf["#CHROM"] == chrname], chrname, nsupp_without_vaf_comparison)
        list_varids_filtout.extend(list_filtout_on_chr)
    
    print("Total filtered out SVs: ", len(list_varids_filtout))
    
    save_vcf_without_specific_varids(path_vcf, path_save, list_varids_filtout)
    

def read_vcf_as_table(path_vcf):
    n_header = check_vcf_header_len(path_vcf)
    table_vcf = pd.read_csv(path_vcf, sep = '\t', skiprows = n_header)
    return table_vcf

def check_vcf_header_len(path_vcf):
    n_header = 0

    read_mode = 'rb' if path_vcf.endswith(".gz") else 'r'
    with open(path_vcf, read_mode) as fr:
        for line in fr:
            if hasattr(line, "decode"):
                line = line.decode()
            if line.startswith("#CHROM"):
                break
            else:
                n_header += 1
    return n_header

def check_variant_overlap_on_single_chromosome(table_vcf_chr, chrname, nsupp_without_vaf_comparison):
    table_vcf_chr_sorted = table_vcf_chr.sort_values(by = "POS")
    
    table_vcf_chr_sorted_del = table_vcf_chr_sorted[table_vcf_chr_sorted["SVType"] == "DEL"]
    table_vcf_chr_sorted_ins = table_vcf_chr_sorted[table_vcf_chr_sorted["SVType"] == "INS"]
    
    list_del_varids_filtin = filter_deletion_svs_with_better_support(table_vcf_chr_sorted_del, nsupp_without_vaf_comparison)
    list_del_varids_filtout = list(filter(lambda varid: varid not in list_del_varids_filtin, table_vcf_chr_sorted_del["ID"]))
    
    list_ins_varids_filtout_overlap = filter_insertion_svs_on_same_position(table_vcf_chr_sorted_ins, nsupp_without_vaf_comparison)
    
    table_vcf_chr_sorted_del_filt = table_vcf_chr_sorted_del[table_vcf_chr_sorted_del["ID"].isin(list_del_varids_filtin)]
    table_vcf_chr_sorted_ins_filt = table_vcf_chr_sorted_ins[~table_vcf_chr_sorted_ins["ID"].isin(list_ins_varids_filtout_overlap)]
    
    list_varids_filtout_indel_overlap = filter_indels_overlap(table_vcf_chr_sorted_del_filt, table_vcf_chr_sorted_ins_filt, nsupp_without_vaf_comparison)
    
    print(chrname, "DEL overlap: ", len(list_del_varids_filtout), "INS overlap: ", len(list_ins_varids_filtout_overlap), "INDEL overlap: ", len(list_varids_filtout_indel_overlap))
    
    list_tot_filtout = list_del_varids_filtout+list_ins_varids_filtout_overlap+list_varids_filtout_indel_overlap
    return list_tot_filtout
    
    
def filter_deletion_svs_with_better_support(table_vcf_chr_sorted_del, nsupp_without_vaf_comparison):
    list_del_varids_filtin = list()
    
    pos_bef = None
    varid_bef = None
    dv_bef = None
    dr_bef = None
    svlen_bef = None
    
    for _, row in table_vcf_chr_sorted_del.iterrows():
        if pos_bef == None:
            pass
        
        elif (pos_bef - svlen_bef ) < row["POS"]:
            list_del_varids_filtin.append(varid_bef)
        
        else:
            if row["DV"] < nsupp_without_vaf_comparison and dv_bef >= nsupp_without_vaf_comparison:
                # print(f"Pass {row['ID']}: DV {row['DV']} vs {dv_bef}")
                continue
            elif row["DV"] >= nsupp_without_vaf_comparison and dv_bef < nsupp_without_vaf_comparison:
                # print(f"Pass {varid_bef}: DV {dv_bef} vs {row['DV']}")
                pass
            else:
                vaf_bef = dv_bef / (dv_bef+dr_bef)
                vaf_curr = row["DV"] / (row["DV"] + row["DR"])
                
                if vaf_bef > vaf_curr:
                    # print(f"Pass {row['ID']}: DV {vaf_curr} vs {vaf_bef}")
                    continue
                elif vaf_bef < vaf_curr:
                    pass
                else:
                    if dv_bef > row["DV"]:
                        continue
                    else:
                        pass
                
        pos_bef = row["POS"]
        varid_bef = row["ID"]
        dv_bef = row["DV"]
        dr_bef = row["DR"]
        svlen_bef = row["SVLEN"]
    
    list_del_varids_filtin.append(varid_bef)
    return list_del_varids_filtin
    
def filter_insertion_svs_on_same_position(table_vcf_chr_sorted_ins, nsupp_without_vaf_comparison):
    count_pos = Counter(table_vcf_chr_sorted_ins["POS"])
    
    list_pos_multiple = list(filter(lambda pos: count_pos[pos] > 1, count_pos.keys()))
    
    list_ins_varids_filtout = list()
    for pos in list_pos_multiple:
        table_vcf_chr_sorted_ins_pos = table_vcf_chr_sorted_ins[table_vcf_chr_sorted_ins["POS"] == pos]
        
        table_vcf_chr_sorted_ins_pos_overdv = table_vcf_chr_sorted_ins_pos[table_vcf_chr_sorted_ins_pos["DV"] >= nsupp_without_vaf_comparison]
        
        if table_vcf_chr_sorted_ins_pos_overdv.shape[0] == 0:
            table_vcf_chr_sorted_ins_pos_overdv = table_vcf_chr_sorted_ins_pos.copy()
            
        table_vcf_chr_sorted_ins_pos_overdv["VAF"] = table_vcf_chr_sorted_ins_pos_overdv["DV"] / (table_vcf_chr_sorted_ins_pos_overdv["DV"]+table_vcf_chr_sorted_ins_pos_overdv["DR"])
        
        max_vaf = table_vcf_chr_sorted_ins_pos_overdv["VAF"].max()
        table_vcf_chr_sorted_ins_pos_overdv_maxvaf = table_vcf_chr_sorted_ins_pos_overdv[table_vcf_chr_sorted_ins_pos_overdv["VAF"] == max_vaf]
        
        table_vcf_chr_sorted_ins_pos_overdv_maxvaf_sorted = table_vcf_chr_sorted_ins_pos_overdv_maxvaf.sort_values(by = "DV", ascending = False)
        
        varid_keep = table_vcf_chr_sorted_ins_pos_overdv_maxvaf_sorted["ID"].to_list()[0]
        
        varids_remove = list(filter(lambda varid: varid != varid_keep, table_vcf_chr_sorted_ins_pos["ID"]))
        list_ins_varids_filtout.extend(varids_remove)
    return list_ins_varids_filtout

def filter_indels_overlap(table_vcf_chr_sorted_del_filt, table_vcf_chr_sorted_ins_filt, nsupp_without_vaf_comparison):
    list_varids_filtout = list()
    
    table_vcf_chr_sorted_del_filt_final = table_vcf_chr_sorted_del_filt.copy()
    table_vcf_chr_sorted_del_filt_final["POSEND"] = table_vcf_chr_sorted_del_filt_final["POS"] - table_vcf_chr_sorted_del_filt_final["SVLEN"]
    
    for _, row_ins in table_vcf_chr_sorted_ins_filt.iterrows():
        for _, row_del in table_vcf_chr_sorted_del_filt_final.iterrows():
            if check_is_value_in_range([row_del["POS"], row_del["POSEND"]], row_ins["POS"]):
                dv_ins = row_ins["DV"]
                dr_ins = row_ins["DR"]
                vaf_ins = dv_ins / (dv_ins + dr_ins)
                
                dv_del = row_del["DV"]
                dr_del = row_del["DR"]
                vaf_del = dv_del / (dv_del + dr_del)
                
                if dv_ins < nsupp_without_vaf_comparison and dv_del >= nsupp_without_vaf_comparison:
                    list_varids_filtout.append(row_ins["ID"])
                    break
                elif dv_ins >= nsupp_without_vaf_comparison and dv_del < nsupp_without_vaf_comparison:
                    list_varids_filtout.append(row_del["ID"])
                    break
                else:
                    if vaf_ins < vaf_del:
                        list_varids_filtout.append(row_ins["ID"])
                    elif vaf_ins > vaf_del:
                        list_varids_filtout.append(row_del["ID"])
                    else:
                        if dv_ins <= dv_del:
                            list_varids_filtout.append(row_ins["ID"])
                        else:
                            list_varids_filtout.append(row_del["ID"])
                            
        table_vcf_chr_sorted_del_filt_final = table_vcf_chr_sorted_del_filt_final[~table_vcf_chr_sorted_del_filt_final["ID"].isin(list_varids_filtout)]
    return list_varids_filtout
        

def check_is_value_in_range(ranges, value):
    return value >= ranges[0] and value <= ranges[1]
    
def save_vcf_without_specific_varids(path_vcf_input, path_vcf_save, list_varids_exclude):
    with open(path_vcf_save, 'w') as fw, open(path_vcf_input, 'r') as fr:
        for line in fr:
            if line.startswith('#'):
                fw.write(line)
            else:
                varid = line.split('\t')[2]
                if varid in list_varids_exclude:
                    continue
                else:
                    fw.write(line)

#%%
if __name__ == "__main__":
    hap = 2
    
    path_vcf = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round4_SV_Call/SV_Call_Merged/Merged_Call/SV_Call.KOREF_PolishedASM_Round4_SV.hap{hap}.hifi_ont.sametrio.filtered.iris.final_evaluated.indel_only.additional_support_of_unclassified_reads_on_depth_5.manual_curation.vcf"
    path_vcf_save = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Variant_Call/Round4_SV_Call/SV_Call_Merged/Merged_Call/SV_Call.KOREF_PolishedASM_Round4_SV.hap{hap}.hifi_ont.sametrio.filtered.iris.final_evaluated.indel_only.additional_support_of_unclassified_reads_on_depth_5.manual_curation.filter_overlap.vcf"
    
    run(
        path_vcf, path_vcf_save, nsupp_without_vaf_comparison=10
    )
# %%
