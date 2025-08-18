#%%
import argparse, os


def run(list_vcfs, path_save):
    with open(path_save, 'w') as fw:
        for ind, path_vcf in enumerate(list_vcfs):
            with open(path_vcf, 'r') as fr:
                for line in fr:
                    if line.startswith('#'):
                        if ind == 0:
                            fw.write(line)
                        else:
                            pass
                    else:
                        fw.write(line)

if __name__ == "__main__":
    argument_parser = argparse.ArgumentParser()

    argument_parser.add_argument("--vcfs", nargs = '+', required = True)
    argument_parser.add_argument("--out")

    args = argument_parser.parse_args()

    print(args)
    run(
        args.vcfs,
        args.out
    )