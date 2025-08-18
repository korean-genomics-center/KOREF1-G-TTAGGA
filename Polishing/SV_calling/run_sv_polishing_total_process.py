#%%
import argparse, sys, os, time, json
from pathlib import Path


sys.path.append(str(Path(__file__).parents[0]))

from run_sniffles import run as run_sniffles
from run_sniffles_with_vcf import run as run_sniffles_on_vcf
from filter_and_refine_svcall_by_iris import run as run_iris
from merge_svcalls_between_hifi_and_ont import run as run_jasmine_hifiont
from merge_svcalls_with_opposite_triobin import run as run_jasmine_difftrio
from run_classify_merged_vcf_with_homocall import run as run_classify_homocall
from run_correct_header_from_genotype_based_sniffles_call import run as run_correct_sniffles_header
from run_filter_homoref_only_from_merged_difftrio_sv import run as run_filter_homoref_from_difftrio
from run_filter_vcf_by_varid_with_prefix import run as run_filter_varid_with_prefix
from run_minimap2_between_assembly_to_assembly import run as run_minimap2_between_asms
from run_paftools_for_assembly_wise_variant_call import run as run_paftools_for_asm_sv
from run_filter_svcall_from_paftools_vcf import run as run_filter_sv_from_paftools
from run_filter_vcf_by_vaf_on_multiple_samples import run as run_filter_vcf_by_vaf_on_multiple_samples
from run_merge_vcf_lines import run as run_merge_vcf_lines
from run_merge_variants_of_jasmine_output_with_rnames import run as run_merge_variants_from_jasmine


THREAD_SNIFFLES2 = 48
THREAD_IRIS = 48
THREAD_JASMINE = 32
THREAD_MINIMAP2 = 48
THREAD_PAFTOOLS = 32

argument_parser = argparse.ArgumentParser()
argument_parser.add_argument("--ref_hap1")
argument_parser.add_argument("--ref_hap2")
argument_parser.add_argument("--bam_hap1_hifi_sametrio")
argument_parser.add_argument("--bam_hap1_hifi_difftrio")
argument_parser.add_argument("--bam_hap1_ont_sametrio")
argument_parser.add_argument("--bam_hap1_ont_difftrio")
argument_parser.add_argument("--bam_hap2_hifi_sametrio")
argument_parser.add_argument("--bam_hap2_hifi_difftrio")
argument_parser.add_argument("--bam_hap2_ont_sametrio")
argument_parser.add_argument("--bam_hap2_ont_difftrio")

argument_parser.add_argument("--dir_output")
argument_parser.add_argument("--name_output")
argument_parser.add_argument("--dir_logging")

argument_parser.add_argument("--json_args")

argument_parser.add_argument("--overwrite", action = "store_true", default = False)

args = argument_parser.parse_args()

if args.json_args != None:
    with open(args.json_args, 'r') as fr:
        dict_args = json.load(fr)
else:
    dict_args = vars(args)


list_job_queue = list()
currtime = int(time.time())

# Run SV calling (Sniffles2)
dir_log_sv = os.path.join(dict_args["dir_logging"], "SV_Calling", "log")
dir_script_sv = os.path.join(dict_args["dir_logging"], "SV_Calling", "sh")
dir_tmp_sv = os.path.join(dict_args["dir_logging"], "SV_Calling", "tmp")

for hap in ["hap1", "hap2"]:
    for tgs in ["hifi", "ont"]:
        for trio in ["sametrio", "difftrio"]:
            if not args.overwrite and os.path.exists(os.path.join(dict_args["dir_output"], "SV_Calling", f"SV_Call.{dict_args['name_output']}.{hap}.{tgs}.{trio}.vcf")):
                continue
            run_sniffles(
                dict_args[f"bam_{hap}_{tgs}_{trio}"],
                dict_args[f"ref_{hap}"],
                os.path.join(dict_args["dir_output"], "SV_Calling", f"SV_Call.{dict_args['name_output']}.{hap}.{tgs}.{trio}.vcf"),
                dir_tmp_sv, dir_script_sv, dir_log_sv, 
                f"SV_Calling.{dict_args['name_output']}.{hap}.{tgs}.{trio}.{currtime}",
                THREAD_SNIFFLES2, jobs_wait = list_job_queue
            )
        
for hap in ["hap1", "hap2"]:
    for tgs in ["hifi", "ont"]:
        for trio in ["sametrio", "difftrio"]:
            list_job_queue.append(f"SV_Calling.{dict_args['name_output']}.{hap}.{tgs}.{trio}.{currtime}")


# Run SV Filtering (filterSV, Iris)
dir_log_iris = os.path.join(dict_args["dir_logging"], "IRIS_SV_Filtering", "log")
dir_script_iris = os.path.join(dict_args["dir_logging"], "IRIS_SV_Filtering", "sh")
dir_tmp_iris = os.path.join(dict_args["dir_logging"], "IRIS_SV_Filtering", "tmp")

for hap in ["hap1", "hap2"]:
    for tgs in ["hifi", "ont"]:
        for trio in ["sametrio", "difftrio"]:
            if not args.overwrite and os.path.exists(os.path.join(dict_args["dir_output"], "SV_Calling", f"SV_Call.{dict_args['name_output']}.{hap}.{tgs}.{trio}.filtered.iris.vcf")):
                continue
            run_iris(
                dict_args[f"ref_{hap}"], 
                os.path.join(dict_args["dir_output"], "SV_Calling", f"SV_Call.{dict_args['name_output']}.{hap}.{tgs}.{trio}.vcf"),
                dict_args[f"bam_{hap}_{tgs}_{trio}"],
                os.path.join(dict_args["dir_output"], "SV_Calling", f"SV_Call.{dict_args['name_output']}.{hap}.{tgs}.{trio}"),
                tgs == "hifi",
                dir_tmp_iris, dir_script_iris, dir_log_iris,
                f"IRIS_SV_Filtering.{dict_args['name_output']}.{hap}.{tgs}.{trio}.{currtime}",
                THREAD_IRIS, jobs_wait=list_job_queue
            )

for tgs in ["hifi", "ont"]:
    for trio in ["sametrio", "difftrio"]:
        list_job_queue.append(f"IRIS_SV_Filtering.{dict_args['name_output']}.{hap}.{tgs}.{trio}.{currtime}")


# Run SV Merging (Jasmine)
dir_log_jasmine = os.path.join(dict_args["dir_logging"], "Jasmine_Merge_SVcall", "log")
dir_script_jasmine = os.path.join(dict_args["dir_logging"], "Jasmine_Merge_SVcall", "sh")
dir_tmp_jasmine = os.path.join(dict_args["dir_logging"], "Jasmine_Merge_SVcall", "tmp")

for hap in ["hap1", "hap2"]:
    if not args.overwrite and os.path.exists(os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.vcf")):
        continue
    run_jasmine_hifiont(
        [os.path.join(dict_args["dir_output"], "SV_Calling", f"SV_Call.{dict_args['name_output']}.{hap}.hifi.sametrio.filtered.iris.vcf"),
         os.path.join(dict_args["dir_output"], "SV_Calling", f"SV_Call.{dict_args['name_output']}.{hap}.ont.sametrio.filtered.iris.vcf")],
        os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris"),
        dir_tmp_jasmine, dir_script_jasmine, dir_log_jasmine,
        f"Jasmine_Merge_SVcall.{dict_args['name_output']}.{hap}.merge_hifiont.sametrio.{currtime}",
        THREAD_JASMINE, jobs_wait = list_job_queue
    )

for hap in ["hap1", "hap2"]:
    list_job_queue.append(f"Jasmine_Merge_SVcall.{dict_args['name_output']}.{hap}.merge_hifiont.sametrio.{currtime}")


# Filter Homo-ALT SVs (Custom)
dir_log_classify_homocall = os.path.join(dict_args["dir_logging"], "Classify_Homo_ALT", "log")
dir_script_classify_homocall = os.path.join(dict_args["dir_logging"], "Classify_Homo_ALT", "sh")
dir_tmp_classify_homocall = os.path.join(dict_args["dir_logging"], "Classify_Homo_ALT", "tmp")

for hap in ["hap1", "hap2"]:
    if not args.overwrite and os.path.exists(os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.hetalt_only.vcf")):
        continue
    run_classify_homocall(
        os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.vcf"),
        os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split"),
        dir_tmp_classify_homocall, dir_script_classify_homocall, dir_log_classify_homocall,
        f"Split_HomoALT_SV.{dict_args['name_output']}.{hap}.merged_hifiont.sametrio.{currtime}",
        1, jobs_wait = list_job_queue
    )

for hap in ["hap1", "hap2"]:
    list_job_queue.append(f"Split_HomoALT_SV.{dict_args['name_output']}.{hap}.merged_hifiont.sametrio.{currtime}")


# Call SVs on Het-ALT regions from Different Trio-bin
for hap in ["hap1", "hap2"]:
    for tgs in ["hifi", "ont"]:
        if not args.overwrite and os.path.exists(os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.hetalt_only.call_difftrio.{tgs}.vcf")):
            continue
        run_sniffles_on_vcf(
            dict_args[f"bam_{hap}_{tgs}_difftrio"],
            dict_args[f"ref_{hap}"],
            os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.hetalt_only.vcf"),
            os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.hetalt_only.call_difftrio.{tgs}.vcf"),
            dir_tmp_sv,
            dir_script_sv,
            dir_log_sv,
            f"SV_Calling.On_HetAlt_Variants.{dict_args['name_output']}.{hap}.{tgs}.difftrio.{currtime}",
            THREAD_SNIFFLES2,
            jobs_wait = list_job_queue
        )

for hap in ["hap1", "hap2"]:
    for tgs in ["hifi", "ont"]:
        list_job_queue.append(f"SV_Calling.On_HetAlt_Variants.{dict_args['name_output']}.{hap}.{tgs}.difftrio.{currtime}")


# Correct Genotype-based Sniffles2 Header
dir_log_correct_header = os.path.join(dict_args["dir_logging"], "Correct_Sniffles_Header", "log")
dir_script_correct_header = os.path.join(dict_args["dir_logging"], "Correct_Sniffles_Header", "sh")
dir_tmp_correct_header = os.path.join(dict_args["dir_logging"], "Correct_Sniffles_Header", "tmp")
for hap in ["hap1", "hap2"]:
    for tgs in ["hifi", "ont"]:
        if not args.overwrite and os.path.exists(os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.hetalt_only.call_difftrio.{tgs}.header_corr.vcf")):
            continue
        run_correct_sniffles_header(
            os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.hetalt_only.call_difftrio.{tgs}.vcf"),
            os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.hetalt_only.call_difftrio.{tgs}.header_corr.vcf"),
            dir_tmp_correct_header,
            dir_script_correct_header,
            dir_log_correct_header,
            f"Correct_Header.On_Sniffles2_Genotype_VCF.{dict_args['name_output']}.{hap}.{tgs}.difftrio.{currtime}",
            1, jobs_wait = list_job_queue
        )

for hap in ["hap1", "hap2"]:
    for tgs in ["hifi", "ont"]:
        list_job_queue.append(f"Correct_Header.On_Sniffles2_Genotype_VCF.{dict_args['name_output']}.{hap}.{tgs}.difftrio.{currtime}")


# Merge Difftrio SV Calls
for hap in ["hap1", "hap2"]:
    if not args.overwrite and os.path.exists(os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.hetalt_only.call_difftrio.total_hifi_ont.merged.vcf")):
        continue
    run_jasmine_difftrio(
        [os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.hetalt_only.call_difftrio.hifi.header_corr.vcf"),
        os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.hetalt_only.call_difftrio.ont.header_corr.vcf"),
        os.path.join(dict_args["dir_output"], "SV_Calling", f"SV_Call.{dict_args['name_output']}.{hap}.hifi.difftrio.filtered.iris.vcf"),
        os.path.join(dict_args["dir_output"], "SV_Calling", f"SV_Call.{dict_args['name_output']}.{hap}.ont.difftrio.filtered.iris.vcf")],
        os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.hetalt_only.call_difftrio.total_hifi_ont"),
        dir_tmp_jasmine, dir_script_jasmine, dir_log_jasmine,
        f"Jasmine_Merge_SVcall.{dict_args['name_output']}.{hap}.merge_opposite_total_triobin_SV_calls.{currtime}",
        THREAD_JASMINE, jobs_wait = list_job_queue, spec_reads = 0
    )

for hap in ["hap1", "hap2"]:
    list_job_queue.append(f"Jasmine_Merge_SVcall.{dict_args['name_output']}.{hap}.merge_opposite_total_triobin_SV_calls.{currtime}")


# Filter HomoRef Variant IDs from Difftrio SV calls
dir_log_filter_homoref_from_difftrio = os.path.join(dict_args["dir_logging"], "Filter_HomoRef_from_Difftrio", "log")
dir_script_filter_homoref_from_difftrio = os.path.join(dict_args["dir_logging"], "Filter_HomoRef_from_Difftrio", "sh")
dir_tmp_filter_homoref_from_difftrio = os.path.join(dict_args["dir_logging"], "Filter_HomoRef_from_Difftrio", "tmp")
for hap in ["hap1", "hap2"]:
    if not args.overwrite and os.path.exists(os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.hetalt_only.call_difftrio.total_hifi_ont.merged.varnames_homoref.txt")):
        continue
    run_filter_homoref_from_difftrio(
        os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.hetalt_only.call_difftrio.total_hifi_ont.merged.vcf"),
        os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.hetalt_only.call_difftrio.total_hifi_ont.merged.varnames_homoref.txt"),
        [0, 1],
        [0, 1, 2, 3],
        dir_tmp_filter_homoref_from_difftrio, dir_script_filter_homoref_from_difftrio, dir_log_filter_homoref_from_difftrio,
        f"Filter_HomoRef_from_Difftrio.{dict_args['name_output']}.{hap}.opposite_total_triobin_SV_calls.{currtime}",
        1, jobs_wait = list_job_queue
    )

for hap in ["hap1", "hap2"]:
    list_job_queue.append(f"Filter_HomoRef_from_Difftrio.{dict_args['name_output']}.{hap}.opposite_total_triobin_SV_calls.{currtime}")


# Filter HomoRef SVs from Difftrio on Het-ALT SVs
dir_log_filter_vcf_with_varids = os.path.join(dict_args["dir_logging"], "Filter_VCF_with_Prefix_Varid", "log")
dir_script_filter_vcf_with_varids = os.path.join(dict_args["dir_logging"], "Filter_VCF_with_Prefix_Varid", "sh")
dir_tmp_filter_vcf_with_varids = os.path.join(dict_args["dir_logging"], "Filter_VCF_with_Prefix_Varid", "tmp")
for hap in ["hap1", "hap2"]:
    if not args.overwrite and os.path.exists(os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.hetalt_only.filt_difftrio_homoref.vcf")):
        continue
    run_filter_varid_with_prefix(
        os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.hetalt_only.vcf"),
        os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.hetalt_only.filt_difftrio_homoref.vcf"),
        os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.hetalt_only.call_difftrio.total_hifi_ont.merged.varnames_homoref.txt"),
        ["0_"],
        dir_tmp_filter_vcf_with_varids, dir_script_filter_vcf_with_varids, dir_log_filter_vcf_with_varids,
        f"Filter_VCF_with_Prefix_Varid.{dict_args['name_output']}.{hap}.difftrio_homoref.{currtime}",
        1, jobs_wait = list_job_queue
    )

for hap in ["hap1", "hap2"]:
    list_job_queue.append(f"Filter_VCF_with_Prefix_Varid.{dict_args['name_output']}.{hap}.difftrio_homoref.{currtime}")


# Align Between Assemblies
dict_hap_exchange = {
    "hap1" : "hap2",
    "hap2" : "hap1"
}
for hap in ["hap1", "hap2"]:
    if not args.overwrite and os.path.exists(os.path.join(dict_args["dir_output"], "ASM_Comparison", f"Align_Minimap2.Ref_{dict_hap_exchange[hap]}.Qry_{hap}.paf")):
        continue
    run_minimap2_between_asms(
        dict_args[f"ref_{hap}"],
        dict_args[f"ref_{dict_hap_exchange[hap]}"],
        os.path.join(dict_args["dir_output"], "ASM_Comparison", f"Align_Minimap2.Ref_{dict_hap_exchange[hap]}.Qry_{hap}.paf"),
        os.path.join(dict_args["dir_logging"], "Minimap2_ASM_Align", "tmp"),
        os.path.join(dict_args["dir_logging"], "Minimap2_ASM_Align", "sh"),
        os.path.join(dict_args["dir_logging"], "Minimap2_ASM_Align", "log"),
        f"Align_Minimap2.{dict_args['name_output']}.Ref_{dict_hap_exchange[hap]}.Qry_{hap}.{currtime}",
        THREAD_MINIMAP2, jobs_wait=list()
    )

for hap in ["hap1", "hap2"]:
    list_job_queue.append(f"Align_Minimap2.{dict_args['name_output']}.Ref_{dict_hap_exchange[hap]}.Qry_{hap}.{currtime}")


# Call SV between ASM
for hap in ["hap1", "hap2"]:
    if not args.overwrite and os.path.exists(os.path.join(dict_args["dir_output"], "ASM_Comparison", f"SV_Call.Paftools.{dict_args['name_output']}.Ref_{dict_hap_exchange[hap]}.Qry_{hap}.svcall.vcf")):
        continue
    run_paftools_for_asm_sv(
        os.path.join(dict_args["dir_output"], "ASM_Comparison", f"Align_Minimap2.Ref_{dict_hap_exchange[hap]}.Qry_{hap}.paf"),
        dict_args[f"ref_{dict_hap_exchange[hap]}"],
        os.path.join(dict_args["dir_output"], "ASM_Comparison", f"SV_Call.Paftools.{dict_args['name_output']}.Ref_{dict_hap_exchange[hap]}.Qry_{hap}.vcf"),
        os.path.join(dict_args["dir_logging"], "Paftools_SV", "tmp"),
        os.path.join(dict_args["dir_logging"], "Paftools_SV", "sh"),
        os.path.join(dict_args["dir_logging"], "Paftools_SV", "log"),
        f"Paftools_SV.{dict_args['name_output']}.Ref_{dict_hap_exchange[hap]}.Qry_{hap}.{currtime}",
        THREAD_PAFTOOLS, jobs_wait=[f"Align_Minimap2.{dict_args['name_output']}.Ref_{dict_hap_exchange[hap]}.Qry_{hap}.{currtime}"]
    )

    run_filter_sv_from_paftools(
        os.path.join(dict_args["dir_output"], "ASM_Comparison", f"SV_Call.Paftools.{dict_args['name_output']}.Ref_{dict_hap_exchange[hap]}.Qry_{hap}.vcf"),
        os.path.join(dict_args["dir_output"], "ASM_Comparison", f"SV_Call.Paftools.{dict_args['name_output']}.Ref_{dict_hap_exchange[hap]}.Qry_{hap}.svcall.vcf"),
        os.path.join(dict_args["dir_logging"], "Paftools_SV", "tmp"),
        os.path.join(dict_args["dir_logging"], "Paftools_SV", "sh"),
        os.path.join(dict_args["dir_logging"], "Paftools_SV", "log"),
        f"Filter_SV_from_Paftools.{dict_args['name_output']}.Ref_{dict_hap_exchange[hap]}.Qry_{hap}.{currtime}",
        THREAD_PAFTOOLS, jobs_wait=[f"Align_Minimap2.{dict_args['name_output']}.Ref_{dict_hap_exchange[hap]}.Qry_{hap}.{currtime}", f"Paftools_SV.{dict_args['name_output']}.Ref_{dict_hap_exchange[hap]}.Qry_{hap}.{currtime}"]

    )

for hap in ["hap1", "hap2"]:
    list_job_queue.append(f"Paftools_SV.{dict_args['name_output']}.Ref_{dict_hap_exchange[hap]}.Qry_{hap}.{currtime}")
    list_job_queue.append(f"Filter_SV_from_Paftools.{dict_args['name_output']}.Ref_{dict_hap_exchange[hap]}.Qry_{hap}.{currtime}")


# Merge Between ASM Calls
for hap in ["hap1", "hap2"]:
    if not args.overwrite and os.path.exists(os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.hetalt_only.filt_difftrio_homoref.sv_between_asm.merged.vcf")):
        continue
    run_jasmine_difftrio(
        [os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.hetalt_only.filt_difftrio_homoref.vcf"),
        os.path.join(dict_args["dir_output"], "ASM_Comparison", f"SV_Call.Paftools.{dict_args['name_output']}.Ref_{hap}.Qry_{dict_hap_exchange[hap]}.svcall.vcf")],
        os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.hetalt_only.filt_difftrio_homoref.sv_between_asm"),
        dir_tmp_jasmine, dir_script_jasmine, dir_log_jasmine,
        f"Jasmine_Merge_SVcall.{dict_args['name_output']}.{hap}.merge_asm_wise_sv.{currtime}",
        THREAD_JASMINE, jobs_wait = list_job_queue, spec_reads = 0
    )

for hap in ["hap1", "hap2"]:
    list_job_queue.append(f"Jasmine_Merge_SVcall.{dict_args['name_output']}.{hap}.merge_asm_wise_sv.{currtime}")


# Filter NA Variant IDs from ASM-wise SV calls
for hap in ["hap1", "hap2"]:
    if not args.overwrite and os.path.exists(os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.hetalt_only.filt_difftrio_homoref.sv_between_asm.merged.varnames_na.txt")):
        continue
    run_filter_homoref_from_difftrio(
        os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.hetalt_only.filt_difftrio_homoref.sv_between_asm.merged.vcf"),
        os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.hetalt_only.filt_difftrio_homoref.sv_between_asm.merged.varnames_na.txt"),
        [],
        [2],
        dir_tmp_filter_homoref_from_difftrio, dir_script_filter_homoref_from_difftrio, dir_log_filter_homoref_from_difftrio,
        f"Filter_HomoRef_from_Difftrio.{dict_args['name_output']}.{hap}.asm_wise_sv.{currtime}",
        1, jobs_wait = list_job_queue
    )

for hap in ["hap1", "hap2"]:
    list_job_queue.append(f"Filter_HomoRef_from_Difftrio.{dict_args['name_output']}.{hap}.asm_wise_sv.{currtime}")


# Filter Missing SVs from ASM-SV on Het-ALT SVs
dir_log_filter_vcf_with_varids = os.path.join(dict_args["dir_logging"], "Filter_VCF_with_Prefix_Varid", "log")
dir_script_filter_vcf_with_varids = os.path.join(dict_args["dir_logging"], "Filter_VCF_with_Prefix_Varid", "sh")
dir_tmp_filter_vcf_with_varids = os.path.join(dict_args["dir_logging"], "Filter_VCF_with_Prefix_Varid", "tmp")
for hap in ["hap1", "hap2"]:
    if not args.overwrite and os.path.exists(os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.hetalt_only.filt_difftrio_homoref.filt_no_asmsv.vcf")):
        continue
    run_filter_varid_with_prefix(
        os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.hetalt_only.filt_difftrio_homoref.vcf"),
        os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.hetalt_only.filt_difftrio_homoref.filt_no_asmsv.vcf"),
        os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.hetalt_only.filt_difftrio_homoref.sv_between_asm.merged.varnames_na.txt"),
        ["0_"],
        dir_tmp_filter_vcf_with_varids, dir_script_filter_vcf_with_varids, dir_log_filter_vcf_with_varids,
        f"Filter_VCF_with_Prefix_Varid.{dict_args['name_output']}.{hap}.no_asmsv.{currtime}",
        1, jobs_wait = list_job_queue
    )

for hap in ["hap1", "hap2"]:
    list_job_queue.append(f"Filter_VCF_with_Prefix_Varid.{dict_args['name_output']}.{hap}.no_asmsv.{currtime}")


# Filter VAF > 0.5, Min_Support > 5 from Final Het-ALT SVs
for hap in ["hap1", "hap2"]:
    if not args.overwrite and os.path.exists(os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.hetalt_only.filt_difftrio_homoref.filt_no_asmsv.vaf_0.5.minsupp_5.vcf")):
        continue
    run_filter_vcf_by_vaf_on_multiple_samples(
        os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.hetalt_only.filt_difftrio_homoref.filt_no_asmsv.vcf"),
        os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.hetalt_only.filt_difftrio_homoref.filt_no_asmsv.vaf_0.5.minsupp_5.vcf"),
        0.5, 
        5,
        os.path.join(dict_args["dir_logging"], "Filter_VAF_Multiple", "tmp"),
        os.path.join(dict_args["dir_logging"], "Filter_VAF_Multiple", "sh"),
        os.path.join(dict_args["dir_logging"], "Filter_VAF_Multiple", "log"),
        f"Filter_VAF_Multiple.{dict_args['name_output']}.{hap}.vaf_0.5.minsupp_5.{currtime}",
        1, jobs_wait=list_job_queue
    )

for hap in ["hap1", "hap2"]:
    list_job_queue.append(f"Filter_VAF_Multiple.{dict_args['name_output']}.{hap}.vaf_0.5.minsupp_5.{currtime}")


# Merge HomALT calls and Filtered HetALT calls
for hap in ["hap1", "hap2"]:
    if not args.overwrite and os.path.exists(os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.homoalt_include.hetalt_only.filt_difftrio_homoref.filt_no_asmsv.vaf_0.5.minsupp_5.vcf")):
        continue
    run_merge_vcf_lines(
        [os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.homoalt_include.vcf"),
        os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.hetalt_only.filt_difftrio_homoref.filt_no_asmsv.vaf_0.5.minsupp_5.vcf")],
        os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.homoalt_include.hetalt_only.filt_difftrio_homoref.filt_no_asmsv.vaf_0.5.minsupp_5.vcf"),
        os.path.join(dict_args["dir_logging"], "Merge_VCF_lines", "tmp"),
        os.path.join(dict_args["dir_logging"], "Merge_VCF_lines", "sh"),
        os.path.join(dict_args["dir_logging"], "Merge_VCF_lines", "log"),
        f"Merge_VCF_lines.{dict_args['name_output']}.{hap}.homalt_filthetalt.{currtime}",
        1, jobs_wait=list_job_queue
    )

for hap in ["hap1", "hap2"]:
    list_job_queue.append(f"Merge_VCF_lines.{dict_args['name_output']}.{hap}.homalt_filthetalt.{currtime}")


# Merge two (HiFi, ONT) SV calls into single SV call
for hap in ["hap1", "hap2"]:
    if not args.overwrite and os.path.exists(os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.homoalt_include.hetalt_only.filt_difftrio_homoref.filt_no_asmsv.vaf_0.5.minsupp_5.call_merged.vcf")):
        continue
    run_merge_variants_from_jasmine(
        os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.homoalt_include.hetalt_only.filt_difftrio_homoref.filt_no_asmsv.vaf_0.5.minsupp_5.vcf"),
        os.path.join(dict_args["dir_output"], "SV_Calling", "Merged_SV", f"SV_Call.{dict_args['name_output']}.{hap}.hifi_ont.sametrio.filtered.iris.merged.split.homoalt_include.hetalt_only.filt_difftrio_homoref.filt_no_asmsv.vaf_0.5.minsupp_5.call_merged.vcf"),
        [os.path.join(dict_args["dir_output"], "SV_Calling", f"SV_Call.{dict_args['name_output']}.{hap}.hifi.sametrio.filtered.iris.vcf"),
         os.path.join(dict_args["dir_output"], "SV_Calling", f"SV_Call.{dict_args['name_output']}.{hap}.ont.sametrio.filtered.iris.vcf")],
        ["RNAMES_HIFI", "RNAMES_ONT"],
        os.path.join(dict_args["dir_logging"], "Merge_Jasmine_Calls", "tmp"),
        os.path.join(dict_args["dir_logging"], "Merge_Jasmine_Calls", "sh"),
        os.path.join(dict_args["dir_logging"], "Merge_Jasmine_Calls", "log"),
        f"Merge_Jasmine_Calls.{dict_args['name_output']}.{hap}.{currtime}",
        1, jobs_wait=list_job_queue
    )

for hap in ["hap1", "hap2"]:
    list_job_queue.append(f"Merge_Jasmine_Calls.{dict_args['name_output']}.{hap}.{currtime}")


