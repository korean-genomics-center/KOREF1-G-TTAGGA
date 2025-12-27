import os, time, json
import subprocess

PATH_SAMTOOLS = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/Samtools/samtools-1.17/samtools-1.17/bin/samtools"
PATH_BEDTOOLS = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/Bedtools/bedtools"


def run_cmd(cmd):
    subprocess.run(cmd, shell = True)

def get_bedtools_intersect_cmd(bed_a, bed_b, path_save, mode = "intersect"):
    cmd = f"{PATH_BEDTOOLS} {mode} -a {bed_a} -b {bed_b} > {path_save}"
    return cmd

def run(path_wholegenome_bed, list_cntr_bed, path_cntrnoct_bed, path_sd_bed, path_sex_bed, dir_out):    
    os.makedirs(dir_out, exist_ok=True)

    dict_annot = dict()

    path_save_biased_annot = os.path.join(dir_out, "biased_annotation_names.txt")
    with open(path_save_biased_annot, 'w') as fw:
        for path_cntr_bed in list_cntr_bed:
            name_save = f"CNTR_{os.path.basename(path_cntr_bed).split('.')[-2]}"
            filename_bed = f"{name_save}.bed"
            path_save = os.path.join(dir_out, filename_bed)
            cmd_intersect = get_bedtools_intersect_cmd(path_cntr_bed, path_wholegenome_bed, path_save)
            run_cmd(cmd_intersect)
            fw.write(name_save + '\n')
            dict_annot[name_save] = f"/annotation/{filename_bed}"

    run_cmd(f"cp {path_wholegenome_bed} {dir_out}/WHOLE_GENOME.bed")
    dict_annot["WHOLE_GENOME"] = "/annotation/WHOLE_GENOME.bed"

    run_cmd(get_bedtools_intersect_cmd(path_cntrnoct_bed, path_wholegenome_bed, os.path.join(dir_out, "CNTR.bed")))
    dict_annot["CNTR"] = "/annotation/CNTR.bed"
    run_cmd(get_bedtools_intersect_cmd(path_wholegenome_bed, path_cntrnoct_bed, os.path.join(dir_out, "Non_CNTR.bed"), mode = "subtract"))
    dict_annot["Non_CNTR"] = "/annotation/Non_CNTR.bed"

    run_cmd(get_bedtools_intersect_cmd(path_sd_bed, path_wholegenome_bed, os.path.join(dir_out, "SD.bed")))
    dict_annot["SD"] = "/annotation/SD.bed"
    run_cmd(get_bedtools_intersect_cmd(path_wholegenome_bed, path_sd_bed, os.path.join(dir_out, "Non_SD.bed"), mode = "subtract"))
    dict_annot["Non_SD"] = "/annotation/Non_SD.bed"

    run_cmd(get_bedtools_intersect_cmd(path_sex_bed, path_wholegenome_bed, os.path.join(dir_out, "SEX.bed")))
    dict_annot["SEX"] = "/annotation/SEX.bed"
    run_cmd(get_bedtools_intersect_cmd(path_wholegenome_bed, path_sex_bed, os.path.join(dir_out, "AUTOSOME.bed"), mode = "subtract"))
    dict_annot["AUTOSOME"] = "/annotation/AUTOSOME.bed"

    run_cmd(
        f"""cat {os.path.join(dir_out, "CNTR.bed")} {os.path.join(dir_out, "SD.bed")} | {PATH_BEDTOOLS} sort -i - | {PATH_BEDTOOLS} merge -i - > {os.path.join(dir_out, "CNTR_or_SD.bed")}"""
    )
    dict_annot["CNTR_or_SD"] = "/annotation/CNTR_or_SD.bed"
    run_cmd(get_bedtools_intersect_cmd(os.path.join(dir_out, "Non_CNTR.bed"), os.path.join(dir_out, "Non_SD.bed"), os.path.join(dir_out, "Non_CNTR_and_Non_SD.bed")))
    dict_annot["Non_CNTR_and_Non_SD"] = "/annotation/Non_CNTR_and_Non_SD.bed"

    run_cmd(get_bedtools_intersect_cmd(os.path.join(dir_out, "AUTOSOME.bed"), os.path.join(dir_out, "Non_CNTR.bed"), os.path.join(dir_out, "AUTOSOME_Non_CNTR.bed")))
    dict_annot["AUTOSOME_Non_CNTR"] = "/annotation/AUTOSOME_Non_CNTR.bed"
    run_cmd(get_bedtools_intersect_cmd(os.path.join(dir_out, "AUTOSOME.bed"), os.path.join(dir_out, "Non_SD.bed"), os.path.join(dir_out, "AUTOSOME_Non_SD.bed")))
    dict_annot["AUTOSOME_Non_SD"] = "/annotation/AUTOSOME_Non_SD.bed"

    run_cmd(get_bedtools_intersect_cmd(os.path.join(dir_out, "AUTOSOME.bed"), os.path.join(dir_out, "Non_CNTR_and_Non_SD.bed"), os.path.join(dir_out, "AUTOSOME_Non_CNTR_and_Non_SD.bed")))
    dict_annot["AUTOSOME_Non_CNTR_and_Non_SD"] = "/annotation/AUTOSOME_Non_CNTR_and_Non_SD.bed"

    with open(os.path.join(dir_out, "annotations_path.json"), 'w') as fw:
        json.dump(dict_annot, fw, indent = 4)


if __name__ == "__main__":
    path_wholegenome_bed = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Flagger_Evaluation/Annotations/Suffix_Attatched_BED/WholeGenome.both.suffixed.bed"
    
    list_cntr_bed = [
        "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Flagger_Evaluation/Annotations/Suffix_Attatched_BED/CenSat.both.suffixed.decomposed.active_hor.bed",
        "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Flagger_Evaluation/Annotations/Suffix_Attatched_BED/CenSat.both.suffixed.decomposed.bsat.bed",
        "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Flagger_Evaluation/Annotations/Suffix_Attatched_BED/CenSat.both.suffixed.decomposed.hsat1A.bed",
        "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Flagger_Evaluation/Annotations/Suffix_Attatched_BED/CenSat.both.suffixed.decomposed.hsat1B.bed",
        "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Flagger_Evaluation/Annotations/Suffix_Attatched_BED/CenSat.both.suffixed.decomposed.hsat2.bed",
        "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Flagger_Evaluation/Annotations/Suffix_Attatched_BED/CenSat.both.suffixed.decomposed.hsat3.bed"
    ]

    path_cntrnoct_bed = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Flagger_Evaluation/Annotations/Suffix_Attatched_BED/CenSat.both.suffixed.no_ct.bed"
    path_sd_bed = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Flagger_Evaluation/Annotations/Suffix_Attatched_BED/SD.both.suffixed.range_merged.bed"
    path_sex_bed = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Flagger_Evaluation/Annotations/Suffix_Attatched_BED/Sex.bed"

    dir_out = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Assembly_Evaluation/Flagger_Evaluation/Annotations/Final_Annotations"
    run(path_wholegenome_bed, list_cntr_bed, path_cntrnoct_bed, path_sd_bed, path_sex_bed, dir_out)