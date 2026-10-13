#%%
import os,sys,math,subprocess

PATH_BEDTOOLS = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/Bedtools/bedtools"


def run_rename(path_bed, path_save, suffix, list_sexchr = ["chrX", "chrY"]):
    with open(path_bed, 'r') as fr, open(path_save, 'w') as fw:
        for line in fr:
            values = line.strip('\n').split('\t')
            if len(values) >= 3:
                chrname = values[0]
                if chrname not in list_sexchr:
                    chrname = f"{chrname}{suffix}"
                values[0] = chrname

                fw.write('\t'.join(values)+'\n')

def run_merge(path_bed1, path_bed2, path_save):
    subprocess.run(f"cat {path_bed1} {path_bed2} > {path_save}", shell = True)

def filter_sexchr(path_bed, path_save, list_sexchr = ["chrX", "chrY"]):
    with open(path_bed, 'r') as fr, open(path_save, 'w') as fw:
        for line in fr:
            chrname = line.split()[0]
            if chrname in list_sexchr:
                fw.write(line)

def run_merge_ranges(path_bed, path_save):
    subprocess.run(f"cat {path_bed} | {PATH_BEDTOOLS} sort -i - | {PATH_BEDTOOLS} merge -i - > {path_save}", shell = True)

#%%
if __name__ == "__main__":
    dict_path_annotations_format = {
        "CenSat" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Repeat_Annotation/Centromere_Annotation/AlphaAnnotation/Hap{hap}/FinalASM.KOREF_S1.Trio_Specific.hap{hap}.Reordered.fasta.rm.cenSat.bed",
        "SD" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Flagger_Evaluation/Annotations/Raw_BED/SD.hap{hap}.err_l10.bed",
        "WholeGenome" : "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Flagger_Evaluation/Annotations/Raw_BED/whole_genome.hap{hap}.bed"
    }

    dir_output = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Flagger_Evaluation/Annotations/Suffix_Attatched_BED"


    for name, path_format in dict_path_annotations_format.items():
        for hap in [1, 2]:
            path_bed = path_format.format(hap = hap)
            output_name = f"{name}.hap{hap}.suffixed.bed"
            path_save = os.path.join(dir_output, output_name)

            suffix = f"_Hap{hap}"

            run_rename(path_bed, path_save, suffix)
        
        run_merge(
            os.path.join(dir_output, f"{name}.hap1.suffixed.bed"),
            os.path.join(dir_output, f"{name}.hap2.suffixed.bed"),
            os.path.join(dir_output, f"{name}.both.suffixed.bed")
        )

        if name == "WholeGenome":
            filter_sexchr(os.path.join(dir_output, f"{name}.both.suffixed.bed"), os.path.join(dir_output, f"Sex.bed"))
        
        if name == "SD":
            run_merge_ranges(os.path.join(dir_output, f"{name}.both.suffixed.bed"), os.path.join(dir_output, f"{name}.both.suffixed.range_merged.bed"))