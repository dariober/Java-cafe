#!/bin/bash

# These tests are far from comprehensive...

stvExe=~/Tritume/SamTextViewer.jar

## MEMO: less -R -S to keep colour
cd ~/svn_git/Java-cafe/trunk/SamTextViewer/test_data/

java -Xmx500m -jar $stvExe -h &&

echo "CAN PRINT FROM ONE POSITION"
java -Xmx500m -jar $stvExe ds051.actb.bam -ni &&

java -Xmx500m -jar $stvExe \
	ds051.actb.bam ds051.actb.bam \
	-fa chr7.fa \
	-r chr7:5566860 -m 10 -ni &&
	
echo "HANDLE NO READS IN INPUT"
java -Xmx500m -jar $stvExe ds051.actb.bam -r chr7:5566860 -m 10 -f 16 -F 16 -ni &&

java -Xmx500m -jar $stvExe ds051.actb.bam -r chr7:5566860 -m 10 -f 16 -ni &&

echo "BED FILES"
java -Xmx500m -jar $stvExe refSeq.hg19.short.bed -ni &&
java -Xmx500m -jar $stvExe refSeq.hg19.short.sort.bed.gz -ni &&

java -Xmx500m -jar $stvExe test.bedGraph.gz -ni &&

echo -e "\n\nDONE\n\n"

exit
