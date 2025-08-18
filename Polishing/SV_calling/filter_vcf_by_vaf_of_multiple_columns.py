import os, argparse


def run(path_vcf, path_save, min_vaf, min_supp, must_all):
    os.makedirs(os.path.dirname(path_save), exist_ok = True)

    with open(path_vcf, 'r') as fr, open(path_save, 'w') as fw:
        for line in fr:
            if line.startswith('#'):
                fw.write(line)
            else:
                list_values = line.strip('\n').split('\t')

                data_format = list_values[8].split(':')
                list_data = list_values[9:]

                list_is_data_satisfy = list()
                for data in list_data:
                    dict_data = dict(zip(data_format, data.split(':')))
                    
                    gt = dict_data["GT"]
                    rs = dict_data["DR"]
                    vs = dict_data["DV"]

                    if gt == "0/0" or gt == "./.":
                        continue

                    rs, vs = int(rs), int(vs)
                    rat = (vs/(rs+vs))

                    if rat >= min_vaf and vs >= min_supp:
                        list_is_data_satisfy.append(True)
                    else:
                        list_is_data_satisfy.append(False)

                if must_all:
                    is_satisfy = sum(list_is_data_satisfy) == len(list_is_data_satisfy)
                else:
                    is_satisfy = sum(list_is_data_satisfy) > 0

                if is_satisfy:
                    fw.write(line)

if __name__ == "__main__":
    argument_parser = argparse.ArgumentParser()

    argument_parser.add_argument("--vcf")
    argument_parser.add_argument("--out")
    argument_parser.add_argument("--min_vaf", type = float, default = 0)
    argument_parser.add_argument("--min_supp", type = int, default = 0)
    argument_parser.add_argument("--must_all", action = "store_true", default = False)

    args = argument_parser.parse_args()
    print(args)

    run(
        args.vcf,
        args.out,
        args.min_vaf,
        args.min_supp,
        args.must_all
    )