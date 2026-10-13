1. run_sv_polishing_total_process.py
Run for both "Trio-Specific only" & "Trio-Specific + Ambiguous" for Same-trio Read based SV calling

2. run_mosdepth.py
Run for bam file with "Trio-Specific only"

3. run_check_depth_of_vcf_variants.py
Check alignment depth of SVs from "Trio-Specific + Ambiguous" dataset on "Trio-Specific only" bam file

4. add_svcalls_only_available_from_unclassified_reads.py

5. filter_sv_with_manual_curation.py
Check large SV calls (> 500bp) on IGV and verify

6. solve_sv_overlap.py