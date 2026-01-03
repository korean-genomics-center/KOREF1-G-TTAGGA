import os, time
import subprocess

PATH_BEDTOOLS = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/Bedtools/bedtools"
PATH_SCRIPT_FILTGFF_BY_GENEONLYGFF = os.path.join(os.path.dirname(__file__), "filter_gff_by_geneonly_parents.py")
PATH_AGAT_GXF2GXF = "/BiO/Access/yoonsung/anaconda3/envs/liftoff/bin/agat_convert_sp_gxf2gxf.pl"


def get_context_for_qsub_submission(cmd, jobname, stderr_path, stdout_path, nthreads, hostname = list(), hold_jids = list()):
    #Header
    script_header = f"#!bin/sh\n"
    script_header += f"#$ -N {jobname}\n"
    script_header += f"#$ -S /bin/bash\n"
    script_header += f"#$ -V\n"
    if len(hold_jids) > 0:
        script_header += f"#$ -hold_jid {','.join(hold_jids)}\n"
    if len(hostname) > 0:
        script_header += f"#$ -l hostname={'|'.join(hostname)}\n"
    script_header += f"#$ -e {stderr_path}\n"
    script_header += f"#$ -o {stdout_path}\n"
    script_header += f"#$ -pe smp {nthreads}\n"
    script_header +="\n"
    script_header+=f"echo \"=================================================\"\n"
    script_header+=f"echo \"Started on	   : $(date)\"\n"
    script_header+=f"echo \"Running on node  : $(hostname)\"\n"
    script_header+=f"echo \"Current job ID   : $JOB_ID\"\n"
    script_header+=f"echo \"Current job name : $JOB_NAME\"\n"
    script_header+=f"echo \"=================================================\"\n\n"

    script_footer="\n\n"
    script_footer+=f"echo \"=================================================\"\n"
    script_footer+=f"echo \"Finished on	   : $(date)\"\n"
    script_footer+=f"echo \"=================================================\"\n"
    script_footer+=f"echo \"Done.\"\n"

    script = script_header + cmd + script_footer
    return script

def write_shell_script_and_run(path_shell, shell_context):
    with open(path_shell, 'w') as fw:
        fw.write(shell_context)
    subprocess.run(f"qsub {path_shell}", shell = True)

def get_commands_to_get_nonoverlaping_genes(path_gff_base, name_base_gff, path_gff_liftoff, name_liftoff_gff, path_output_prefix):
    cmd = f"""
awk '$3=="gene"' {path_gff_base} > {path_output_prefix}.{name_base_gff}_gene.gff3
awk '$3=="gene"' {path_gff_liftoff} > {path_output_prefix}.{name_liftoff_gff}_gene.gff3
{PATH_BEDTOOLS} intersect -v -a {path_output_prefix}.{name_liftoff_gff}_gene.gff3 -b {path_output_prefix}.{name_base_gff}_gene.gff3 > {path_output_prefix}.{name_liftoff_gff}_gene.nonoverlap.gff3
python {PATH_SCRIPT_FILTGFF_BY_GENEONLYGFF} --gffin {path_gff_liftoff} --gffgeneonly {path_output_prefix}.{name_liftoff_gff}_gene.nonoverlap.gff3 --gffout {path_output_prefix}.{name_liftoff_gff}.nonoverlap.gff3

{PATH_BEDTOOLS} intersect -a {path_output_prefix}.{name_liftoff_gff}.nonoverlap.gff3 -b {path_output_prefix}.{name_base_gff}_gene.gff3 > {path_output_prefix}.{name_liftoff_gff}.evaluate_nonoverlap.txt

cat {path_gff_base} {path_output_prefix}.{name_liftoff_gff}.nonoverlap.gff3 > {path_output_prefix}.{name_base_gff}_{name_liftoff_gff}.gff3
{PATH_AGAT_GXF2GXF} -gff {path_output_prefix}.{name_base_gff}_{name_liftoff_gff}.gff3 -o {path_output_prefix}.{name_base_gff}_{name_liftoff_gff}.sorted.gff3

"""
    return cmd, f"{path_output_prefix}.{name_base_gff}_{name_liftoff_gff}.sorted.gff3"

def run(path_cat_gff, dict_liftoff_gffs, path_output_prefix, dir_tmp, dir_shell, dir_log, jobname, n_jobs, hostname = list(), jobs_wait = list()):    
    path_gff_base = path_cat_gff
    name_gff_base = "CAT"

    cmd = ""
    for name_liftoff, path_gff_liftoff in dict_liftoff_gffs.items():
        cmd_add_nonoverlap, path_result_gff = get_commands_to_get_nonoverlaping_genes(path_gff_base, name_gff_base, path_gff_liftoff, name_liftoff, path_output_prefix)
        cmd += cmd_add_nonoverlap

        path_gff_base = path_result_gff
        name_gff_base = f"{name_gff_base}__{name_liftoff}"

    cmd += f"ln -s {path_result_gff} {path_output_prefix}.{name_gff_base}.liftoff_final.gff3"
        
    os.makedirs(os.path.dirname(path_output_prefix), exist_ok=True)
    os.makedirs(dir_tmp, exist_ok=True)
    os.makedirs(dir_shell, exist_ok=True)
    os.makedirs(dir_log, exist_ok=True)
    
    write_shell_script_and_run(
        f"{dir_shell}/{jobname}.sh",
        get_context_for_qsub_submission(
            cmd,
            jobname,
            f"{dir_log}/{jobname}.stderr",
            f"{dir_log}/{jobname}.stdout",
            n_jobs,
            hostname=hostname,
            hold_jids=jobs_wait
        )
    )



if __name__ == "__main__":
    dir_workspace = "/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Workspace/Yoonsung/KOREF_Assembly/Analysis"
    
    curr_time = int(time.time())

    for hap in [1, 2]:
        path_cat_gff = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Run_CAT/Hap{hap}.Blood_Cellline.CGP_PB/consensus_gene_set/KOREF1_ASM_Hap{hap}.gff3"

        dict_liftoff_gffs = {
            "CHM13_Liftoff":f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/Liftoff/From.chm13v2_CATLiftoff_v2.To.KOREF_ASM_20251027_Hap{hap}/liftoff_annotation.KOREF_ASM_Hap{hap}.from_chm13v2_CATLiftoff_v2.annotate.gff_polished",
            "GRCh38_Liftoff":f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/Liftoff/From.GRCh38_GENCODEv49.To.KOREF_ASM_20251027_Hap{hap}/liftoff_annotation.KOREF_ASM_Hap{hap}.from_GRCh38_GENCODEv49.annotate.gff3_polished"
        }

        path_output_prefix = f"/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Gene_Annotation/CAT_Annotation/Postprocess_CAT/Hap{hap}.Blood_Cellline.CGP_PB.Add_Liftoff_CHM13_GRCh38/KOREF1_ASM_Hap{hap}"

        for path in dict_liftoff_gffs.values():
            assert os.path.exists(path)

        assert os.path.exists(path_cat_gff)
        run(
            path_cat_gff, dict_liftoff_gffs, path_output_prefix, 
            f"{dir_workspace}/CAT_Annotation/Add_Nonoverlap_Liftoff/tmp", 
            f"{dir_workspace}/CAT_Annotation/Add_Nonoverlap_Liftoff/sh", 
            f"{dir_workspace}/CAT_Annotation/Add_Nonoverlap_Liftoff/log", 
            f"Add_Nonoverlap_Liftoff.KOREF1_FinalASM.CAT.CHM13_GRCh38.hap{hap}.{curr_time}",
            1,
            jobs_wait=[]
        )