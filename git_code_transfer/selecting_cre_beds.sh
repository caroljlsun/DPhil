#!/bin/bash

#---------------------------------------------
# Procedure: select CREs from encode files in usgc
# Project: Doctoral work
# User: CZJS
# Date: August_2023
#---------------------------------------------

#echo "Running!"


#bigBedToBed http://hgdownload.soe.ucsc.edu/gbdb/hg38/encode3/ccre/encodeCcreCombined.bb -chrom=chr21 -start=0 -end=100000000 stdout

#bigBedToBed http://hgdownload.soe.ucsc.edu/gbdb/hg38/encode3/ccre/encodeCcreCombined.bb -chrom=chr10 -start=122098811 -end=122099629 stdout

while read line; do
  chromosome=$(echo $line | cut -d ' ' -f 1)
  start=$(echo $line | cut -d ' ' -f 2)
  stop=$(echo $line | cut -d ' ' -f 3)
  
  
  bigBedToBed http://hgdownload.soe.ucsc.edu/gbdb/hg38/encode3/ccre/encodeCcreCombined.bb -chrom=$chromosome -start=$start -end=$stop stdout
  
  #echo $chromosome $start $stop
done

#echo "Stopping!"

#grep '\-1' cre_intersections.bed