import os, time
import subprocess

PATH_SCRIPT_FILTGFF_BY_GENEONLYGFF = os.path.join(os.path.dirname(os.path.dirname(__file__)), "filter_gff_by_geneonly_parents.py")


def get_commands_to_get_nonoverlaping_genes(path_gff, path_output_prefix):
    cmd = f"""
grep "unknown_likely_coding" {path_gff} | awk '$3=="gene" {{print $0}}' > {path_output_prefix}.unknown_gene.gff3 &&
python {PATH_SCRIPT_FILTGFF_BY_GENEONLYGFF} --gffin {path_gff} --gffgeneonly {path_output_prefix}.unknown_gene.gff3 --gffout {path_output_prefix}.unknown_full.gff3 &&
awk '$3=="transcript" {{print $0}}' {path_output_prefix}.unknown_full.gff3 > {path_output_prefix}.unknown_transcript.gff3
"""
    return cmd

def run(path_gff, path_output_prefix):    
    cmd = get_commands_to_get_nonoverlaping_genes(path_gff, path_output_prefix)
        
    os.makedirs(os.path.dirname(path_output_prefix), exist_ok=True)
    
    subprocess.run(cmd, shell = True)



if __name__ == "__main__":

    for hap in [1, 2]:
        path_gff = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Run_CAT/Hap{hap}.Blood_Cellline.CGP_PB/consensus_gene_set/KOREF1_ASM_Hap{hap}.gff3"

        path_output_prefix = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/Process_Unknown/Hap{hap}.Blood_Cellline.CGP_PB/KOREF1_ASM_Hap{hap}"

        run(path_gff, path_output_prefix)