1. run_make_transcript_fasta.py
run gffread to make cDNA, Amino acid fasta files
(Whole gff & Unknown-only gff)

2. run_diamond_on_protein_fasta.py
Input: Unknown-only protein CDS fasta
-> filter evalue <= 1e-5

3. check_diamond_interproscan_result.py

4. make_salmon_index.py

5. 