import os


def read_gff_as_dict(path_gff):
    dict_gff = dict()

    with open(path_gff, 'r') as fr:
        for line in fr:
            if line.startswith('#'): continue

            


def get_attrib_values_from_gff_line(line):
    attrib = line.strip('\n').split('\t')[8]
    values = attrib.split(';')

    dict_attrib = dict()
    for val in values:
        if len(val.split('=')) != 2: continue
        key, key_val = val.split('=')
        dict_attrib[key] = key_val

    return dict_attrib