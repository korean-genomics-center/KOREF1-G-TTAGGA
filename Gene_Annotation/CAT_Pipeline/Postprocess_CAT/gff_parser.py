import os


def read_gff_as_dict(path_gff):
    dict_gff = {
        "metainfo" : list(),
        "data" : list()
    }

    general_structure = ["seqid", "source", "feat_type", "start", "end", "score", "strand", "phase"]

    with open(path_gff, 'r') as fr:
        for line in fr:
            if line.startswith('#'): 
                dict_gff["metainfo"].append(line)
            else:
                values = line.strip('\n').split('\t')
                assert len(values) == 9, f"Not proper structure: {line}"

                dict_line_info = dict(zip(general_structure, values[:8]))
                dict_line_info["attrib"] = get_attrib_values_from_gff_line(values[8])

                dict_gff["data"].append(dict_line_info)
    return dict_gff

def get_attrib_values_from_gff_line(attrib_value):
    # attrib_value = line.strip('\n').split('\t')[8]
    values = attrib_value.split(';')

    dict_attrib = dict()
    for val in values:
        if len(val.split('=')) != 2: 
            print(f"Skip: {val}")
            continue
        key, key_val = val.split('=')
        dict_attrib[key] = key_val

    return dict_attrib

def save_dict_gff_into_file(dict_gff, path_save):
    general_structure = ["seqid", "source", "feat_type", "start", "end", "score", "strand", "phase"]

    with open(path_save, 'w') as fw:
        fw.write(''.join(dict_gff["metainfo"]))

        for dict_value in dict_gff["data"]:
            list_values_write = list(map(dict_value.__getitem__, general_structure))
            
            list_attrib_values = list()
            for att_key, att_val in dict_value["attrib"].items():
                list_attrib_values.append(f"{att_key}={att_val}")
            list_values_write.append(';'.join(list_attrib_values))

            fw.write('\t'.join(list_values_write)+'\n')