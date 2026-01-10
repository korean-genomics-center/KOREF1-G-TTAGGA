#%%
import argparse, os
from itertools import chain


def run(path_gff, path_gff_gene_only, path_save):
    os.makedirs(os.path.dirname(path_save), exist_ok = True)
    list_ids_filtin = get_ids_from_gff(path_gff_gene_only)

    set_ids_filtin = set(list_ids_filtin)
    assert len(list_ids_filtin) == len(set_ids_filtin)

    with open(path_gff, 'r') as fr, open(path_save, 'w') as fw:
        for line in fr:
            if line.startswith('#'):
                # fw.write(line)
                pass
            else:
                dict_attrib = get_attrib_values_from_gff_line(line)
                featid = dict_attrib.get("ID", '')
                parid = dict_attrib.get("Parent", '')
                
                if featid in set_ids_filtin:
                    fw.write(line)
                elif parid in set_ids_filtin:
                    set_ids_filtin.add(featid)
                    fw.write(line)

def get_ids_from_gff(path_gff, key_find = "ID"):
    list_ids = list()

    with open(path_gff, 'r') as fr:
        for line in fr:
            id_parsed = False
            if line.startswith('#'):
                continue
            else:
                attribs = line.strip('\n').split('\t')[8]
                for attrib in attribs.split(';'):
                    values = attrib.split('=')
                    if len(values) == 2:
                        key, val = values
                        if key == key_find:
                            assert not id_parsed, line
                            list_ids.append(val)
                            id_parsed = True
            assert id_parsed, line
    return list_ids

def get_attibs_filtin(lineid):
    geneline_id = f"ID={lineid}"
    others_child_id = f"gene={lineid}"
    return geneline_id, others_child_id

def get_attrib_values_from_gff_line(line):
    attrib = line.strip('\n').split('\t')[8]
    values = attrib.split(';')

    dict_attrib = dict()
    for val in values:
        if len(val.split('=')) != 2: continue
        key, key_val = val.split('=')
        dict_attrib[key] = key_val

    return dict_attrib

#%%
if __name__ == "__main__":
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument("--gffin")
    argument_parser.add_argument("--gffgeneonly")
    argument_parser.add_argument("--gffout")
    # argument_parser.add_argument("--gene_find_attrib_key", default = "ID")

    args = argument_parser.parse_args()

    run(
        args.gffin, args.gffgeneonly, args.gffout
    )
