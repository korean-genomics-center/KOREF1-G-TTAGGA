#%%
import argparse, os


def run(path_vcf, path_save):
    os.makedirs(os.path.dirname(path_save), exist_ok = True)

    with open(path_vcf, 'r') as fr, open(path_save, 'w') as fw:
        for line in fr:
            if line.startswith("#CHROM"):
                break
            else:
                fw.write(line)
        header = line
        list_header_values = header.strip('\n').split('\t')

        first_gt_line = fr.readline()
        list_gt_values = first_gt_line.strip('\n').split('\t')

        assert len(list_gt_values) <= len(list_header_values)

        new_header = '\t'.join(list_header_values[:len(list_gt_values)]) + '\n'

        fw.write(new_header)
        fw.write(first_gt_line)

        for line in fr:
            fw.write(line)


if __name__ == "__main__":
    argument_parser = argparse.ArgumentParser()

    argument_parser.add_argument("--vcf")
    argument_parser.add_argument("--out")

    args = argument_parser.parse_args()

    run(
        args.vcf,
        args.out
    )