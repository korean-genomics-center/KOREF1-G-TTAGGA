#!/usr/bin/env bash

echo "Usage: spectra-cn.sh [-c] <read.meryl> <asm1.fasta> [asm2.fasta] out-prefix"
echo -e "\t-c\t\t: [OPTIONAL] input meryl databases are homopolymer compressed, while asm.fasta isn't"
echo -e "\t<read.meryl>\t: Generated with meryl count from i.e. illumina wgs reads"
echo -e "\t<asm1.fasta>\t: haplotype 1 assembly. gzipped or not"
echo -e "\t[asm2.fasta]\t: haplotype 2 assembly. gzipped or not"
echo -e "\t<out-prefix>: output prefix. Required."
echo -e "\t\tWhen only <asm1.fasta> is given, results will be generated in haploid mode."
echo -e "\t\tWhen <asm2.fasta> is given, results will be generated for each asm1 asm2 haploid assembly and asm1+asm2 diploid assembly."
echo


if [[ $# -lt 3 ]]; then
  echo "No args provided. Exit."
  exit -1
fi

source $MERQURY/util/util.sh

compress=""
if [ "x$1" = "x-c" ]; then
  compress="compress"
  echo
  echo "Compress mode enabled."
  echo "  spectra-cn, spectra-asm plots and qv will be generated."
  echo "  for the error .bed and .wig tracks,"
  echo "  use pre-compressed assemblies with no -c option."
  echo
  shift
fi

read=`link $1`
asm1_fa=`link $2`
name=$4

k=`meryl print $read | head -n 2 | tail -n 1 | awk '{print length($1)}'`
echo "Detected k-mer size $k"
echo

if [ -z $name ]; then
  name=$3
  asm2_fa=""
else
  asm2_fa=`link $3`
fi

if [ -s $name ]; then
  echo "$name already exists. Provide a different name."
  exit -1
fi

has_module=$(check_module)
if [[ $has_module -gt 0 ]]; then
  echo "No modules available.."
else
  module load R
fi

asm1=`echo $asm1_fa | sed 's/.fa$//g' | sed 's/.fasta$//g' | sed 's/.fasta.gz$//g' | sed 's/.fa.gz$//g'`


for asm_fa in $asm1_fa $asm2_fa	# will generate only for asm1_fa if asm2_fa is empty
do
  asm=`echo $asm_fa | sed 's/.fa$//g' | sed 's/.fasta$//g' | sed 's/.fasta.gz$//g' | sed 's/.fa.gz$//g'`

  if [ ! -e $asm.k$k.meryl ]; then
    echo "# Generate meryl db for $asm"
    meryl count k=$k $compress output ${asm}.k$k.meryl $asm_fa
    echo
  fi

  echo "# Copy numbers in k-mers found only in asm"
  meryl difference output $asm.k$k.0.meryl ${asm}.k$k.meryl $read

  echo "# QV statistics"
  ASM_ONLY=`meryl statistics ${asm}.k$k.0.meryl  | head -n4 | tail -n1 | awk '{print $2}'`
  TOTAL=`meryl statistics ${asm}.k$k.meryl  | head -n4 | tail -n1 | awk '{print $2}'`
  ERROR=`echo "$ASM_ONLY $TOTAL" | awk -v k=$k '{print (1-(1-$1/$2)^(1/k))}'`
  QV=`echo "$ASM_ONLY $TOTAL" | awk -v k=$k '{print (-10*log(1-(1-$1/$2)^(1/k))/log(10))}'`
  echo -e "$asm\t$ASM_ONLY\t$TOTAL\t$QV\t$ERROR" > $name.qv
  echo

  echo "# Per seq QV statistics"
  meryl-lookup -existence -sequence $asm_fa -mers $asm.k$k.0.meryl/ | \
  awk -v k=$k '{print $1"\t"$4"\t"$2"\t"(-10*log(1-(1-$4/$2)^(1/k))/log(10))"\t"(1-(1-$4/$2)^(1/k))}' > $name.perseq.qv
  echo
done