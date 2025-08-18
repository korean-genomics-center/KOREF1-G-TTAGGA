#%%
import os, argparse

def run(path_vcf, path_save):
    with open(path_vcf, 'r') as fr, open(path_save, 'w') as fw:
        for line in fr:
            if line.startswith('#'):
                fw.write(line)
            else:
                list_values = line.split('\t')
                
                refcall = list_values[3]
                altcall = list_values[4]
                
                if set(list(refcall)) <= {'A', 'T', 'G', 'C'}:
                    pass
                else:
                    print(f"Skip REF: {refcall} in {line}")
                    continue
                
                if set(list(altcall)) <= {'A', 'T', 'G', 'C'}:
                    pass
                else:
                    print(f"Skip ALT: {altcall} in {line}")
                    continue
                
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