1. run_pepperdeepvariant.py
HiFi: pm0a, {same-trio}0a, {same-trio}, {diff-trio}
ONT: {same-trio}0a, {same-trio}, {diff-trio}

2. filter_pass_variant.py
HiFi: pm0a, {same-trio}0a, {same-trio}, {diff-trio}
ONT: {same-trio}0a, {same-trio}, {diff-trio}

3. filter_trio_binned_vcfs.py
HiFi: Passed {same-trio}0a, Passed {same-trio}, Passed {diff-trio}
ONT: Passed {same-trio}0a, Passed {same-trio}, Passed {diff-trio}

4. split_and_filter_vcf_by_biallelic_and_multiallelic.py
HiFi: Passed pm0a

5. filter_ont_trio_binned_homoalt_vcfs.py
ONT: Passed {same-trio}, Passed {same-trio}0a

6. run_mosdepth.py
Run for bam file with "Trio-Specific only"

7. run_check_depth_of_vcf_variants.py
ONT: Filtered Biallelic Homozygous Passed {same-trio}, Filtered Biallelic Homozygous Passed {same-trio}0a, modepth {same-trio}
ONT: Passed {same-trio}, Passed {same-trio}0a, modepth {same-trio}
HiFi: Passed {same-trio}, Passed {same-trio}0a, modepth {same-trio}

8. filter_snvcalls_with_unclassified_by_depth_on_specific.py
depth: 5
ONT: Filtered Biallelic Homozygous Passed {same-trio}0a, depth_on_variants {same-trio}
ONT: Passed {same-trio}0a, depth_on_variants {same-trio}
HiFi: Passed {same-trio}0a, depth_on_variants {same-trio}

7. evaluate_biallelic_het_snv.py
HiFi: Filtered Biallelic Heterozygous Passed pm0a
HiFi: Filtered Passed {same-trio}, Filtered Passed Unseen {same-trio}0a, Filtered Passed {diff-trio}, gVCF {diff-trio}
ONT: Filtered Passed {same-trio}, Filtered Passed Unseen {same-trio}0a, Filtered Passed {diff-trio}, gVCF {diff-trio}

8. evaluate_multiallelic_het_snv.py
HiFi: Filtered Multiallelic Passed pm0a
HiFi: Filtered Passed {same-trio}, Filtered Passed Unseen {same-trio}0a, Filtered Passed {diff-trio}
ONT: Filtered Passed {same-trio}, Filtered Passed Unseen {same-trio}0a, Filtered Passed {diff-trio}

9. merge_vcfs_for_polishing.py
HiFi: Filtered Biallelic Homozygous Passed pm0a
HiFi: Evaluated Filtered Biallelic Heterozygous Passed pm0a
HiFi: Evaluated Filtered Multiallelic Passed pm0a
ONT: Filtered Biallelic Homozygous Passed {same-trio}
ONT: Filtered Biallelic Homozygous Passed Unseen {same-trio}0a