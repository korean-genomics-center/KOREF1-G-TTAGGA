/BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Resources/Tools/Samtools/samtools-1.17/samtools-1.17/samtools view -F 260 /BiO/Research/Project1/KOREF_PersonalMultiomicsReference/Results/KOREF_Assembly/Polishing/Read_Alignment/Align.Winnowmap2/Align_Gapclosed_ASM/Primary_Alignment/Winnowmap_Align.Primary.Ref.KOREF_ASM.hap2.QRY.ONT_UL.Filter_L10kb_QV15.triobin_split.triobin_m.sorted.bam | \
awk 'BEGIN {total_identity=0; count=0} {
    for (i=12; i<=NF; ++i) {
        if ($i ~ /^NM:i:/) {
            nm=substr($i, 6);
            rlen=length($10);
            if (rlen > 0) {
                identity = (rlen - nm) / rlen;
                total_identity += identity;
                 count++;
             }
             break;
         }
     }
 } END {
     if (count > 0) {
         printf "Processed %d alignments.\n", count;
         printf "Average Identity: %.4f\n", total_identity / count;
     } else {
         print "No alignments with NM tag found or processed.";
     }
 }'