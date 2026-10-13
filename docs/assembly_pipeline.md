# KOREF1-G-TTAGGA assembly pipeline

Detailed notes on the scripts used to build the **KOREF1 telomere-to-telomere diploid assembly** (KOREF_S1, trio-phased), organised by pipeline stage.

The code was written and run by Yoonsung Kwon on the KOGIC SGE cluster. The commit history up to 2026-03-29 is preserved from the original repository (`KOGIC_KOREF_AGGA`). The stage layout below was added afterwards by moving files; script contents are unchanged except for `sys.path` lines (see [Notes](#notes)).

## Final assembly lineage

```
01 Input filtering   HiFi (L≥1 kb, QV≥20) · ONT ultra-long (L≥100 kb, QV≥15) · parental WGS → k31 parent-specific hapmers
        │
02 Initial assembly  hifiasm 0.25.0 trio mode (-1/-2 parent-specific yak) + ONT UL      → Trio_Specific hap1 / hap2
        │
03 Scaffolding       tidk telomere split → T2T contigs assigned to chromosomes
        │            → ntLink on non-T2T contigs (trio-binned ONT) → RagTag to CHM13 → chromosome-level scaffolds
        │
04 Gap filling       modified TGS-GapCloser, 9 rounds (manual scaffolding, local / binned assemblies,
        │            unplaced scaffolds, trio-binned HiFi reads)                    → Final.Gapclosed_Assembly
        │
05 Polishing         R1–R4: SV (Sniffles2/IRIS/Jasmine) → SNV (PEPPER-DeepVariant + Merfin)
        │            R5–R7: manual replacement of low-QV windows
        │            R8:    SV (Merfin -polish) → SNV                              → Round8.SV_SNV_Polished_Fasta
        │
06 Final assembly    copy of Round8 hap1/hap2 · Combined assembly by per-chromosome vote · Diploid FASTA
```

The final haplotypes (`Final_Assembly/FinalASM.KOREF_S1.Trio_Specific.hap{1,2}.fasta`) are byte-identical to `Polishing/Reference_Fasta/Round8.SV_SNV_Polished_Fasta`.

## Repository layout

| Folder | Contents |
|---|---|
| `01_Input_Filtering/` | Read preprocessing and trio k-mer databases |
| `02_Initial_Assembly/` | hifiasm trio assembly, GFA → FASTA, assembly stats |
| `03_Scaffolding/` | Telomere-aware scaffolding without Hi-C (`a`–`e` in run order) |
| `04_Gap_Filling/` | TGS-GapCloser rounds and the donor sequences they used |
| `05_Polishing/` | Read alignment, k-mer DB, SV / SNV / manual polishing |
| `06_Final_Assembly/` | Combined and diploid FASTA generation |
| `common/` | Shared helpers (`fasta_parser.py`, modified Merqury scripts) |
| `Evaluation/` | QC of the final assembly (Merqury, GCI, CRAQ, Flagger, Compleasm, NucFreq, Hi-C, figures) |
| `Annotation/` | Gene annotation (CAT, Liftoff) and repeat annotation (RepeatMasker, BISER, centromere) |
| `Archive/` | Superseded or exploratory branches, kept for the record |

## Stages

Most `run_*.py` files are job generators: they build a command, write an SGE script to `Analysis/<Step>/<Sample>/sh/`, and submit it with `qsub`. Inputs and outputs are set in each script's `__main__` block.

### 01 Input filtering

| Step | Scripts | Tool / key parameters |
|---|---|---|
| ONT adapter trim | `ONT/run_dorado_trim.py` | dorado 0.9.6 `trim --sequencing-kit SQK-ULK114` |
| ONT end trim | `ONT/trim_fastq.py` | 1,500 bp removed from each end |
| HiFi CCS → FASTQ, adapter filter | `HiFi/convert_ccs_bam_to_fastq.py`, `HiFi/filter_hifi_adapter.py` | bam2fastq, HiFiAdapterFilt |
| Length / quality filter | `Read_Filter/filter_read_with_stats.py`, `get_read_stats_as_pickle.py` | HiFi L≥1,000 QV≥20; ONT UL L≥100 kb QV≥15 (L≥50 kb QV≥20 for ntLink, L≥10 kb QV≥15 for polishing) |
| Merge runs | `Merge/run_merge_fqgz.py`, `merge_fqgz_files.py` | |
| Parental WGS trim | `Trio_Hapmer/run_fastp_trimming.py` | fastp `--trim_front 15 --trim_tail 5 --n_base_limit 0`, mean/average quality 20 |
| Parent-specific k-mers | `Trio_Hapmer/run_meryl_kmer_count.py`, `run_hapmer.py`, `run_meryl_print_kmer_count.py`, `generate_kmer_fasta_from_merqury_print.py`, `run_yak_kmer_count_single_fastq.py` | meryl k=31 → Merqury `hapmers.sh` → `yak count -b0 -k31` |
| Read trio-binning | `Trio_Hapmer/run_yak_kmer_count.py`, `split_fastq_by_triobin_sort.py` | `yak count -b37 -k31`, `yak triobin` |

### 02 Initial assembly

`run_hifiasm.py` (hifiasm 0.25.0):
`--dual-scaf --telo-m CCCTAA -1 <paternal>.yak -2 <maternal>.yak --ul <ONT L100kb QV15>` on HiFi L1000 QV20 reads.
Output: `Initial_Assembly/Hifiasm/KOREF_S1.ul_L100kb_QV15.Trio_Specific.default_option`.
Then `convert_gfa_into_fasta.py` and `run_gfastats_for_asseembly_stats.py`.

### 03 Scaffolding (no Hi-C)

| Folder | What it does |
|---|---|
| `a_Telomere_Split/` | tidk search for `CCCTAA`; splits contigs into telomere-both / telomere-single (masked) / no-telomere |
| `b_T2T_Chromosome_Assign/` | minimap2 to CHM13, PAF filtering, MUMmer plots; assigns T2T contigs to chromosomes |
| `c_ntLink/` | ntLink `scaffold gap-fill k=32 w=100 sensitive=True` on non-T2T contigs with trio-binned ONT L50kb reads; linkage adjusted around masked telomeres |
| `d_RagTag_CHM13/` | RagTag `scaffold -r -w` to CHM13 v2.0 (minimap2 2.29); reorders telomere-end scaffolds |
| `e_Chromosome_Level/` | Builds chromosome-level scaffolds (100 N gaps); the `n_pseudo_telorepeat5` version was carried forward |

### 04 Gap filling

Gap closing used a modified TGS-GapCloser: <https://github.com/Yoonsung1203/TGS-GapCloser>. The KOREF1 gap filling ran between 2025-08-18 and 2025-09-02, when the latest commit was `dd8b556` (2025-07-25). Commits from 2025-12-21 onward came after the final assembly.

| Round | Donor / method | Scripts |
|---|---|---|
| R1 | Scaffold split + manual scaffolding from a JSON config | `TGS_GapCloser/run_tgsseqsplit_for_split_scaffold.py`, `Manual_Scaffolding/run_manual_scaffolding_by_config.py` |
| R2 | Local HiFi assemblies (hap1 chr13, hap2 chr15) | `Local_Assembly_Donor/*`, `TGS_GapCloser/run_tgsgapcloser_with_assemblies.py` |
| R3 | Unplaced scaffolds | `TGS_GapCloser/run_tgsgapcloser_with_unplaced_scaffolds.py` |
| R4–R7 | Trio-binned hifiasm / Flye assemblies, contigs filtered to QV≥50 and Hamming ≤0.05 | `Binned_Assembly_Donor/*`, `common/Merqury_modified/filter_contig_by_{qv,hamming}.py`, `run_tgsgapcloser_with_assemblies.py` |
| R8 | Manual scaffolding with custom PAF alignments | `Manual_Scaffolding/*` (incl. `run_tgsgapcloser_with_custom_alignment.py`) |
| R9 | Trio-binned HiFi reads (hap1 `-x lr:hq -w5 -e0`, idy_hard 0.97; hap2 `-x map-hifi`, idy_hard 0.98) | `TGS_GapCloser/run_tgsgapcloser_with_HiFi.py` |

TGS-GapCloser settings shared across rounds: `--ne --tgstype pb --min_idy 0.9 --hard`.

### 05 Polishing

| Folder | Role |
|---|---|
| `Alignment/` | Winnowmap 2.03 (fork) alignment of trio-binned HiFi (`map-pb`) and ONT L10kb (`map-ont`) reads to each round's assembly; primary-only filter, merging, mosdepth |
| `Kmer_DB/` | Meryl k=21 read k-mers and GenomeScope2 fitted histogram for Merfin |
| `SV_calling/` | `run_sv_polishing_total_process.py` (driver, `--json_args config.json`): Sniffles2 → IRIS → Jasmine (HiFi+ONT, opposite trio-bin) → paftools cross-check → filters; then manual curation of large SVs and `run_bcftools_consensus.py` (`-H A`). Step order: `README.txt` |
| `SNV_calling/` | PEPPER-Margin-DeepVariant r0.8 → PASS / trio / bi- & multi-allelic filters → depth checks → merge → Merfin 1.1 `-polish` (k=21) → `bcftools consensus`. Step order: `README.txt` |
| `Manual_Polishing/` | Rounds 5–7: find QV<60 windows, align binned-assembly contigs anchored on both sides, check QV, Hamming and uniqueness, replace the reference segment (extra_range 25 / 100 / 250 bp). Step order: `README.txt` |

Round order: gap-closed → R1–R4 (SV then SNV each round) → R5, R6, R7 (manual) → R8 SV (Merfin `-polish` output kept) → R8 SNV → final.

### 06 Final assembly

- `gather_assembly_informations.py`: per-chromosome vote between haplotypes on length, Merqury QV, GCI, CRAQ AQI and Flagger reliable ratio (ties go to the longer one). Output: `Vote_Chromosomes_for_Combined.tsv`.
- `generate_combined_fasta.py`: builds the Combined assembly from the vote.
- `generate_diploid_fasta.py`: concatenates hap1 + hap2 into the diploid FASTA used for Flagger.

## Steps with no script in this repository

- **Upper-casing the gap-closed FASTA, copying Round8 → `Final_Assembly`, chromosome reordering (`*.Reordered.fasta`):** done by hand on the command line.
- **Mitochondrial genome:** assembled in a separate Korea10K pipeline and copied in.
- **Input filtering and hifiasm runs before 2025-05-21:** these ran before the first commit. The scripts' `__main__` blocks were edited afterwards (for example `run_hifiasm.py` now points at a later `Trio_Hapmer` test). The exact commands that were run are in the cluster's `Analysis/*/sh/` job scripts.
- **External tools:** the modified TGS-GapCloser (linked above), IRIS `filter.py` from CallSV, and Merqury `trio/hapmers.sh`.

## Notes

- **Paths are server-specific.** Input, output, tool and `dir_workspace` paths are hard-coded to `/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/...` and were left as they were run. Some scripts reused from other projects still point to their paths: `Annotation/Repeat_Annotation/Satellite_Annotation/*` (Endangered_Jindo_Dog), `Archive/HiC_Scaffolding/*` and `Archive/Local_Assembly/Assemble_missed_telomere_end/*` (Endangered_Goral).
- **Imports.** Helper modules now live in `common/`. Scripts that used `sys.path.append(Path(__file__).parents[N])` to reach the old repository root now append `<repo>/common` instead. That one line is the only content change made during the reorganisation.
- **Generators that call other scripts by absolute path** still use the original server layout (e.g. `.../Scripts/Merqury_modified/get_qv_only.sh`).
- **Not every file in a lineage folder was used in the final round.** `05_Polishing/SV_calling`, `SNV_calling` and `04_Gap_Filling/Manual_Scaffolding` were kept whole so that sibling imports keep working. They may contain helpers that only earlier rounds used.
