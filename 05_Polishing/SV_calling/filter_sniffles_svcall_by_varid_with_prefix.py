#%%
import argparse, os


def run(path_vcf, path_varid, path_save, prefix_remove):
    os.makedirs(os.path.dirname(path_save), exist_ok = True)

    with open(path_varid, 'r') as fr:
        list_varids = fr.readlines()
    list_varids = list(map(lambda val: val.strip('\n'), list_varids))

    list_varids_without_prefix = list(map(lambda val: remove_prefix_from_value(val, prefix_remove), list_varids))

    set_varids_without_prefix = set(list_varids_without_prefix)

    list_varids_found = list()
    with open(path_vcf, 'r') as fr, open(path_save, 'w') as fw:
        for line in fr:
            if line.startswith('#'):
                fw.write(line)
            else:
                list_values = line.strip('\n').split('\t')
                varid = list_values[2]
                if varid in set_varids_without_prefix:
                    fw.write(line)
                    list_varids_found.append(varid)
    
    print(f"{len(set(list_varids_found))} / {len(set_varids_without_prefix)} Variant IDs found from VCF for filtering")

def remove_multiple_prefix_from_value(value, list_prefix):
    while 1:
        removed = False
        for prefix in list_prefix:
            if value.startswith(prefix):
                value = value[len(prefix):]
                removed = True
        
        if not removed:
            break
    return value

def remove_prefix_from_value(value, list_prefix):
    for prefix in list_prefix:
        if value.startswith(prefix):
            value = value[len(prefix):]
            break
    return value


#%%
if __name__ == "__main__":
    argument_parser = argparse.ArgumentParser()

    argument_parser.add_argument("--vcf")
    argument_parser.add_argument("--varid")
    argument_parser.add_argument("--out")
    argument_parser.add_argument("--prefix_remove", default = list(), nargs = '+')

    args = argument_parser.parse_args()

    run(
        args.vcf,
        args.varid,
        args.out,
        args.prefix_remove
    )