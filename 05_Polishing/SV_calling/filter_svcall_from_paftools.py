import os, argparse


def run(path_vcf, path_save, min_svlen = 30):
    os.makedirs(os.path.dirname(path_save), exist_ok = True)

    with open(path_vcf, 'r') as fr, open(path_save, 'w') as fw:
        for line in fr:
            if line.startswith('#'):
                fw.write(line)
            else:
                list_values = line.strip('\n').split('\t')
                rname = list_values[0]
                refseq = list_values[3]
                altseq = list_values[4]
                info = list_values[7]

                list_info_flags, dict_info = parse_info_data_from_vcf(info)

                if dict_info["QNAME"] != rname:
                    continue
                if abs(len(refseq) - len(altseq)) < min_svlen:
                    continue

                fw.write(line)

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


if __name__  == "__main__":
    argument_parser = argparse.ArgumentParser()

    argument_parser.add_argument("--vcf")
    argument_parser.add_argument("--out")
    argument_parser.add_argument("--minlen", type = int, default = 30)
    
    args = argument_parser.parse_args()

    run(
        args.vcf,
        args.out,
        args.minlen
    )