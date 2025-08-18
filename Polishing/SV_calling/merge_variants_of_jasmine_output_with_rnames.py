#%%
import argparse, os
from copy import deepcopy


FM_GT = "GT"
FM_REF = "DR"
FM_ALT = "DV"
FM_VTYPE = "OT"
INFO_RNAMS = "RNAMES"

def run(
    path_vcf_jasmine,
    list_vcfs_merged,
    info_rnames,
    path_save
):
    list_dict_varid_to_rnames = list(map(lambda path_vcf: get_varid_to_support_reads(path_vcf), list_vcfs_merged))
    with open(path_vcf_jasmine, 'r') as fr, open(path_save, 'w') as fw:
        for line in fr:
            if line.startswith("#CHROM"):
                list_header = line.strip('\n').split('\t')
                list_samples = list_header[9:]
                assert len(list_samples) == len(list_vcfs_merged)

                newsampleid = '__'.join(list_samples)
                fw.write('\t'.join(list_header[:9]+[newsampleid])+'\n')

            elif line.startswith('#'):
                fw.write(line)

            else:
                list_values = line.strip('\n').split('\t')
                list_values_write = deepcopy(list_values[:9])
                infoval = list_values[7]
                list_samplevalues = list_values[9:]
                list_info_flags, dict_info_items = parse_info_data_from_vcf(infoval)

                list_suppvec = list(dict_info_items["SUPP_VEC"])
                list_idlist = dict_info_items["IDLIST"].split(',')

                ind_idlist = 0
                for suppvec, dict_varid_to_rnames, rnames_add in zip(list_suppvec, list_dict_varid_to_rnames, info_rnames):
                    if suppvec == '0':
                        dict_info_items[rnames_add] = ''
                    elif suppvec == '1':
                        svid_find = list_idlist[ind_idlist]
                        ind_idlist += 1
                        rnames_on_svid = dict_varid_to_rnames.get(svid_find)
                        assert rnames_on_svid != None, svid_find
                        dict_info_items[rnames_add] = rnames_on_svid
                    else:
                        raise Exception(f"Unsupported SUPP_VEC: {dict_info_items['SUPP_VEC']}")
                list_values_write[7] = get_info_data_from_values(list_info_flags, dict_info_items)

                list_samplevalues_support = list(map(
                    lambda ind: list_samplevalues[ind],
                    list(filter(lambda val: list_suppvec[val] == '1', range(len(list_suppvec))))    
                ))
                variantcall_rep = merge_variant_calls(list_samplevalues_support, list_values[8])
                list_values_write.append(variantcall_rep)

                fw.write('\t'.join(list_values_write) + '\n')

            

def get_varid_to_support_reads(path_vcf):
    dict_varid_to_rnames = dict()

    with open(path_vcf, 'r') as fr:
        for line in fr:
            if line.startswith('#'):
                continue

            list_values = line.strip('\n').split('\t')
            varid = list_values[2]
            infoval = list_values[7]
            list_flags, dict_info_items = parse_info_data_from_vcf(infoval)

            dict_varid_to_rnames[varid] = dict_info_items.get("RNAMES")
    return dict_varid_to_rnames

def parse_info_data_from_vcf(info_data):
    list_data = info_data.split(';')
    
    list_flags = list()
    dict_item = dict()
    
    for data in list_data:
        if '=' in data:
            key, value = data.split('=')
            assert dict_item.get(key) == None
            dict_item[key] = value
        else:
            list_flags.append(data)
    return list_flags, dict_item

def get_info_data_from_values(list_flags, dict_item):
    list_info_data = list()

    list_info_data.extend(list_flags)
    list_info_data.extend(list(map(lambda keyval: f"{keyval[0]}={keyval[1]}", dict_item.items())))

    return ';'.join(list_info_data)

def merge_variant_calls(list_calls, call_fm):
    list_dict_call = list(map(lambda call: dict(zip(call_fm.split(':'), call.split(':'))), list_calls))


    list_gts = list(map(lambda dict_call: dict_call[FM_GT], list_dict_call))
    list_ref = list(map(lambda dict_call: int(dict_call[FM_REF]), list_dict_call))
    list_alt = list(map(lambda dict_call: int(dict_call[FM_ALT]), list_dict_call))
    list_vartype = list(map(lambda dict_call: dict_call[FM_VTYPE], list_dict_call))

    rep_gt = "./."
    if "1/1" in list_gts:
        rep_gt = "1/1"
    elif "0/1" in list_gts:
        rep_gt = "0/1"
    else:
        raise Exception(f"Unchecked Genotypes (1/1 or 0/1 not included): {list_gts}")
    
    sum_ref = sum(list_ref)
    sum_alt = sum(list_alt)
    assert len(set(list_vartype)) == 1, list_vartype
    rep_vartype = list_vartype[0]

    dict_repcall = {
        FM_GT:rep_gt,
        FM_REF:sum_ref,
        FM_ALT:sum_alt,
        FM_VTYPE:rep_vartype
    }

    repcall = ':'.join(
        list(map(str, 
                list(map(lambda fm: dict_repcall.get(fm, '.'), 
                        call_fm.split(':')
                ))
        ))
    )
    return repcall


if __name__ == "__main__":
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument("--vcf_jasmine")
    argument_parser.add_argument("--vcfs_merged", nargs = '+', default = list())
    argument_parser.add_argument("--rnames_info", nargs = '+', default = list())
    argument_parser.add_argument("--out")

    args = argument_parser.parse_args()
    print(args)
    run(
        args.vcf_jasmine,
        args.vcfs_merged,
        args.rnames_info,
        args.out
    )