# KOREF1-G-TTAGGA: Korean Standard Reference Genome

**KOREF** (KOrean REFerence) is the Korean standard reference genome project run by the Korean Genomics Center (KOGIC) at UNIST.

**KOREF1-G-TTAGGA** is the Telomere-to-Telomere Accurate and Gapless Genome Assembly (TTAGGA) of KOREF1, a Korean male donor (sample KOREF-S1). It is a diploid assembly with fully phased paternal and maternal haplotypes, built from PacBio HiFi, Oxford Nanopore ultra-long and parental short reads. Both haplotypes are gapless and telomere-to-telomere for every chromosome. The assembly is integrated with over a decade of multi-omic data from the same person (genome, transcriptome, epigenome, proteome, microbiome and clinical data). Reference materials from the same donor are also available: a B lymphoblastoid cell line (KCLB #60211) and an iPSC line.

All data are open access without any restriction.

- Portal: **<https://koreanreference.org>**
- Dataset page: <https://koreanreference.org/datasets/KOREF1-G-TTAGGA/>
- Genome browser: <https://koreanreference.org/browser/>
- Reference materials: <https://koreanreference.org/materials/>

This repository holds the code used to build the assembly (see [Code in this repository](#code-in-this-repository)).

## Assembly

| Haplotype | Parent | Sequences | Size | Merqury QV |
|---|---|---|---|---|
| KOREF1_HAP1 | Paternal | chr1–22, chrY | 2.914 Gb | 81.19 |
| KOREF1_HAP2 | Maternal | chr1–22, chrX | 3.029 Gb | 79.03 |

Haplotypes were separated by trio binning. Hi-C and other orthogonal libraries were used to evaluate the assembly, not to build it.

## Downloads

All files are served over HTTPS from `koreanreference.org`. FASTA, GFF3, BED and BEDPE files are bgzip-compressed and come with `.fai`/`.gzi` or `.tbi` indexes in the same directory (append the index extension to the file URL). MD5 and SHA-256 checksums for every file are in the manifest on the [dataset page](https://koreanreference.org/datasets/KOREF1-G-TTAGGA/).

### Assembly

- [KOREF1_HAP1.reference.fa.gz](https://koreanreference.org/data/releases/korefome-p0-reference-v1/assemblies/KOREF1_HAP1/reference/reference.fa.gz): paternal haplotype (720 MB). MD5 `0106c3ef57d5e7529997b50d65757473`
- [KOREF1_HAP2.reference.fa.gz](https://koreanreference.org/data/releases/korefome-p0-reference-v1/assemblies/KOREF1_HAP2/reference/reference.fa.gz): maternal haplotype (751 MB). MD5 `f68cab9c3b7f6f4a3066d89044a84fc2`

### Gene annotation

[KOREF1-G-TTAGGA: Personal Genome Annotation](https://koreanreference.org/datasets/koref1-ttagga-genome-annotation/)

- [KOREF1_HAP1.genes.gff3.gz](https://koreanreference.org/data/releases/korefome-p0-reference-v1/assemblies/KOREF1_HAP1/annotation/genes.gff3.gz)
- [KOREF1_HAP2.genes.gff3.gz](https://koreanreference.org/data/releases/korefome-p0-reference-v1/assemblies/KOREF1_HAP2/annotation/genes.gff3.gz)

Gene models from the Comparative Annotation Toolkit (CAT), supplemented with Liftoff of GENCODE (GRCh38) and the CHM13 annotation.

### Repeat annotation

- RepeatMasker
  - [KOREF1_HAP1.repeatmasker.bed.gz](https://koreanreference.org/data/releases/korefome-p1-annotation-tracks-v2/tracks/KOREF1_HAP1/repeats/repeatmasker.bed.gz)
  - [KOREF1_HAP2.repeatmasker.bed.gz](https://koreanreference.org/data/releases/korefome-p1-annotation-tracks-v2/tracks/KOREF1_HAP2/repeats/repeatmasker.bed.gz)
- Centromeric satellites (CenSat)
  - [KOREF1_HAP1.alpha_censat.bed.gz](https://koreanreference.org/data/releases/korefome-p1-annotation-tracks-v2/tracks/KOREF1_HAP1/centromeres/alpha_censat.bed.gz)
  - [KOREF1_HAP2.alpha_censat.bed.gz](https://koreanreference.org/data/releases/korefome-p1-annotation-tracks-v2/tracks/KOREF1_HAP2/centromeres/alpha_censat.bed.gz)
- Active centromeres (live α-satellite HOR arrays)
  - [KOREF1_HAP1.alpha_active.bed.gz](https://koreanreference.org/data/releases/korefome-p1-annotation-tracks-v2/tracks/KOREF1_HAP1/centromeres/alpha_active.bed.gz)
  - [KOREF1_HAP2.alpha_active.bed.gz](https://koreanreference.org/data/releases/korefome-p1-annotation-tracks-v2/tracks/KOREF1_HAP2/centromeres/alpha_active.bed.gz)

### Segmental duplications (BISER)

`strict` and `broad` call sets. Each set has three files: pairwise alignments (`pairs.bedpe`), duplicated intervals (`endpoints.bed`) and merged coverage (`coverage.bed`).

| | strict | broad |
|---|---|---|
| HAP1 | [pairs](https://koreanreference.org/data/releases/korefome-p1-annotation-tracks-v2/tracks/KOREF1_HAP1/segmental_duplications/biser_strict.pairs.bedpe.gz) · [endpoints](https://koreanreference.org/data/releases/korefome-p1-annotation-tracks-v2/tracks/KOREF1_HAP1/segmental_duplications/biser_strict.endpoints.bed.gz) · [coverage](https://koreanreference.org/data/releases/korefome-p1-annotation-tracks-v2/tracks/KOREF1_HAP1/segmental_duplications/biser_strict.coverage.bed.gz) | [pairs](https://koreanreference.org/data/releases/korefome-p1-annotation-tracks-v2/tracks/KOREF1_HAP1/segmental_duplications/biser_broad.pairs.bedpe.gz) · [endpoints](https://koreanreference.org/data/releases/korefome-p1-annotation-tracks-v2/tracks/KOREF1_HAP1/segmental_duplications/biser_broad.endpoints.bed.gz) · [coverage](https://koreanreference.org/data/releases/korefome-p1-annotation-tracks-v2/tracks/KOREF1_HAP1/segmental_duplications/biser_broad.coverage.bed.gz) |
| HAP2 | [pairs](https://koreanreference.org/data/releases/korefome-p1-annotation-tracks-v2/tracks/KOREF1_HAP2/segmental_duplications/biser_strict.pairs.bedpe.gz) · [endpoints](https://koreanreference.org/data/releases/korefome-p1-annotation-tracks-v2/tracks/KOREF1_HAP2/segmental_duplications/biser_strict.endpoints.bed.gz) · [coverage](https://koreanreference.org/data/releases/korefome-p1-annotation-tracks-v2/tracks/KOREF1_HAP2/segmental_duplications/biser_strict.coverage.bed.gz) | [pairs](https://koreanreference.org/data/releases/korefome-p1-annotation-tracks-v2/tracks/KOREF1_HAP2/segmental_duplications/biser_broad.pairs.bedpe.gz) · [endpoints](https://koreanreference.org/data/releases/korefome-p1-annotation-tracks-v2/tracks/KOREF1_HAP2/segmental_duplications/biser_broad.endpoints.bed.gz) · [coverage](https://koreanreference.org/data/releases/korefome-p1-annotation-tracks-v2/tracks/KOREF1_HAP2/segmental_duplications/biser_broad.coverage.bed.gz) |

### Regulatory and CpG annotation

- CpG islands with shores, shelves and open sea
  - [KOREF1_HAP1.cgi_elements.bed.gz](https://koreanreference.org/data/releases/korefome-p1-annotation-tracks-v2/tracks/KOREF1_HAP1/cpg/cgi_elements.bed.gz)
  - [KOREF1_HAP2.cgi_elements.bed.gz](https://koreanreference.org/data/releases/korefome-p1-annotation-tracks-v2/tracks/KOREF1_HAP2/cpg/cgi_elements.bed.gz)
- ENCODE SCREEN v4 candidate cis-regulatory elements (cCREs) placed on KOREF1
  - [KOREF1_HAP1.ccre_screen_v4.bed.gz](https://koreanreference.org/data/releases/korefome-p1-annotation-tracks-v2/tracks/KOREF1_HAP1/regulatory/ccre_screen_v4.bed.gz)
  - [KOREF1_HAP2.ccre_screen_v4.bed.gz](https://koreanreference.org/data/releases/korefome-p1-annotation-tracks-v2/tracks/KOREF1_HAP2/regulatory/ccre_screen_v4.bed.gz)

### Multi-omic data

Sequencing and other omics data from KOREF-S1 (195 experiments across 6 omics layers) are listed in the [data catalog](https://koreanreference.org/catalog/).

## Notes on downloading files

No login is required. Use `wget --continue` (or `curl -C -`) to resume interrupted transfers, for example:

```bash
wget --continue https://koreanreference.org/data/releases/korefome-p0-reference-v1/assemblies/KOREF1_HAP1/reference/reference.fa.gz \
     -O KOREF1_HAP1.reference.fa.gz
```

The dataset pages also provide a `manifest.tsv` and a generated verified-fetch script.

## Code in this repository

The scripts that produced the assembly, ordered by pipeline stage. Most scripts generate and submit SGE jobs, and paths are hard-coded to the KOGIC cluster. Treat them as a record of what was run, not a portable pipeline. Details are in [docs/assembly_pipeline.md](docs/assembly_pipeline.md).

| Folder | Stage |
|---|---|
| `01_Input_Filtering/` | HiFi, ONT ultra-long and parental short-read filtering; parent-specific k-mers |
| `02_Initial_Assembly/` | hifiasm trio assembly with ONT ultra-long reads |
| `03_Scaffolding/` | Telomere-aware scaffolding (ntLink, RagTag to CHM13) |
| `04_Gap_Filling/` | Gap closing with a [modified TGS-GapCloser](https://github.com/Yoonsung1203/TGS-GapCloser) |
| `05_Polishing/` | SV, SNV and manual polishing rounds |
| `06_Final_Assembly/` | Combined and diploid FASTA |
| `07_Evaluation/` | Final-assembly QC (Merqury, Hi-C, GCI, CRAQ, Compleasm, SyRI, Flagger) |
| `08_Gene_Annotation/` | Iso-Seq, CAT and Liftoff gene annotation |
| `09_Repeat_Annotation/` | RepeatMasker, BISER segmental duplications, centromeres |
| `common/` · `Archive/` | Shared helpers, and superseded branches |

## Data reuse and license

All KOREF data are released as open access without any restriction. If you use KOREF1-G-TTAGGA, please cite the reference publication below.

### KOREF1-G-TTAGGA

- Bhak J, Kwon Y, Shin D-H, An K, Choi H, Park S, Yoon C, Choi Y, et al. [Telomere-to-Telomere Accurate and Gapless Korean Standard Reference Genome](https://doi.org/10.1101/2025.11.17.688257). bioRxiv, 2025. *Preprint, not peer reviewed.*

### Earlier KOREF releases

- Kim HS, Jeon S, Kim Y, Kim C, Bhak J, Bhak J. [KOREF_S1: phased, parental trio-binned Korean reference genome using long reads and Hi-C sequencing methods](https://doi.org/10.1093/gigascience/giac022). GigaScience, 2022.
- Kim HM, Jeon S, Chung O, Jun JH, Kim HS, Blazyte A, Lee HY, Yu Y, et al. [Comparative analysis of 7 short-read sequencing platforms using the Korean Reference Genome: MGI and Illumina sequencing benchmark for whole-genome sequencing](https://doi.org/10.1093/gigascience/giab014). GigaScience, 2021.
- Kim HS, Jeon S, Kim C, Kim YK, Cho YS, Kim J, Blazyte A, Manica A, et al. [Chromosome-scale assembly comparison of the Korean Reference Genome KOREF from PromethION and PacBio with Hi-C mapping information](https://doi.org/10.1093/gigascience/giz125). GigaScience, 2019.
- Cho YS, Kim H, Kim HM, Jho S, Jun J, Lee YJ, Chae KS, Kim CG, et al. [An ethnically relevant consensus Korean reference genome is a step towards personal reference genomes](https://doi.org/10.1038/ncomms13637). Nature Communications, 2016.

## Contact

Jong Hwa Bhak, Korean Genomics Center (KOGIC), Ulsan National Institute of Science and Technology (UNIST): jongbhak@genomics.org
