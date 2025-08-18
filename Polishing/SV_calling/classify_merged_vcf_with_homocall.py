#%%
import os, sys, argparse


def run(path_vcf, path_output_prefix):
    os.makedirs(os.path.dirname(path_output_prefix), exist_ok = True)

    path_output_homocall = f"{path_output_prefix}.homoalt_include.vcf"
    path_output_hetcall = f"{path_output_prefix}.hetalt_only.vcf"

    with open(path_vcf, 'r') as fr:
        with open(path_output_homocall, 'w') as fw_homo, open(path_output_hetcall, 'w') as fw_het:
            for line in fr:
                if line.startswith("#"):
                    fw_homo.write(line)
                    fw_het.write(line)
                else:
                    list_values = line.strip('\n').split('\t')
                    list_sample_values = list_values[9:]

                    list_genotypes = list(map(lambda val: val.split(':')[0], list_sample_values))

                    assert set(list_genotypes) <= set({"./.", "0/0", "0/1", "1/1"}), line.strip('\n')

                    if "1/1" in list_genotypes:
                        fw_homo.write(line)
                    else:
                        fw_het.write(line)


if __name__ == "__main__":
    argument_parser = argparse.ArgumentParser()

    argument_parser.add_argument("--vcf")
    argument_parser.add_argument("--out_prefix")

    args = argument_parser.parse_args()

    run(
        args.vcf,
        args.out_prefix
    )