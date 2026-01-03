#%%
import argparse, os
from itertools import chain


def run(path_gff, path_gff_gene_only, path_save):
    os.makedirs(os.path.dirname(path_save), exist_ok = True)
    list_ids_filtin = get_ids_from_gff(path_gff_gene_only)
    assert len(list_ids_filtin) == len(set(list_ids_filtin))

    list_attribs_filtin = list(chain(*
                                     list(map(get_attibs_filtin, list_ids_filtin))))
    set_attribs_filtin = set(list_attribs_filtin)

    with open(path_gff, 'r') as fr, open(path_save, 'w') as fw:
        for line in fr:
            if line.startswith('#'):
                # fw.write(line)
                pass
            else:
                list_attibs = line.strip('\n').split('\t')[8].split(';')
                if len(set(list_attibs) & set_attribs_filtin) > 0:
                    fw.write(line)
                # attrib_values = get_attrib_values_from_gff_line(line)
                # if len(set(attrib_values) & set_ids_filtin) > 0:
                #     fw.write(line)


def get_ids_from_gff(path_gff):
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
                        if key == "ID":
                            assert not id_parsed, line
                            list_ids.append(val)
                            id_parsed = True
            assert id_parsed, line
    return list_ids

def get_attibs_filtin(lineid):
    geneline_id = f"ID={lineid}"
    others_child_id = f"gene_id={lineid}"
    return geneline_id, others_child_id

def get_attrib_values_from_gff_line(line):
    attrib = line.strip('\n').split('\t')[8]
    values = attrib.split(';')

    list_attrib_values = list(map(lambda val: val.split('=')[1], 
                                  list(filter(lambda val2: len(val2.split('=')) == 2, values))
                        ))
    return list_attrib_values


if __name__ == "__main__":
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument("--gffin")
    argument_parser.add_argument("--gffgeneonly")
    argument_parser.add_argument("--gffout")

    args = argument_parser.parse_args()

    run(
        args.gffin, args.gffgeneonly, args.gffout
    )
