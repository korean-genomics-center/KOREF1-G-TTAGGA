#%%
import argparse, os


def run(path_vcf, path_save, inds_must_homoref, inds_check):
    os.makedirs(os.path.dirname(path_save), exist_ok = True)

    list_varnames_homorefcall = list()
    with open(path_vcf, 'r') as fr:
        for line in fr:
            if line.startswith("#"):
                continue
            else:
                list_values = line.strip('\n').split('\t')

                varname = list_values[2]

                list_sample_values = list_values[9:]
                list_gt_values = list(map(lambda val: val.split(':')[0], list_sample_values))

                list_gt_values_must_homoref = list(map(lambda ind: list_gt_values[ind], inds_must_homoref))
                list_gt_values_check = list(map(lambda ind: list_gt_values[ind], inds_check))

                include_this_var = True

                if len(list_gt_values_must_homoref) > 0:
                    if "0/0" not in list_gt_values_must_homoref:
                        include_this_var = False

                if len(list_gt_values_check) > 0:
                    if set(list_gt_values_check) <= {"0/0", "./."}:
                        pass
                    else:
                        include_this_var = False
                
                if include_this_var:
                    list_varnames_homorefcall.append(varname)
    
    with open(path_save, 'w') as fw:
        for varname in list_varnames_homorefcall:
            fw.write(varname + '\n')


if __name__ == "__main__":
    argument_parser = argparse.ArgumentParser()
    
    argument_parser.add_argument("--vcf")
    argument_parser.add_argument("--out")
    argument_parser.add_argument("--inds_homoref", default = list(), nargs = "+", type = int)
    argument_parser.add_argument("--inds_check_nonalt", default = list(), nargs = "+", type = int)

    args = argument_parser.parse_args()

    run(
        args.vcf,
        args.out,
        args.inds_homoref,
        args.inds_check_nonalt
    )