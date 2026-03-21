# Preprocessing
1. run_make_transcript_fasta.py
run gffread to make cDNA, Amino acid fasta files
(Whole gff & Unknown-only gff)

# Protein Homology
1. run_diamond_on_protein_fasta.py
Input: Unknown-only protein CDS fasta
-> filter evalue <= 1e-5

2. run_interprocan_on_protein_fasta.py

3. check_diamond_interproscan_result.py

# Expression Quantification
1. make_salmon_decoy.py

2. make_salmon_index.py

3. run_salmon_expression.py

4. get_transcripts_not_expressed.py

5. evaluate_transcripts_by_salmon.py

# CDS Length
1. get_short_cds_transcripts.py

# Final Evaluation
1. gather_filtering_results.py

2. draw_upset_plot.py

3. get_final_filterout_transcripts.py

4. filter_out_qcfailed_transcripts.py