#!/bin/bash

## MEMO: less -R -S to keep colour
cd /.../SamTextViewer/test

java -Xmx500m -jar ../SamTextViewer.jar -h

echo "CAN PRINT FROM ONE POSITION"
java -Xmx500m -jar ../SamTextViewer.jar -i test_data/ds051.actb.bam -r chr7 

java -Xmx500m -jar ../SamTextViewer.jar \
	-i test_data/ds051.actb.bam test_data/ds051.actb.bam \
	-fa test_data/chr7.fa \
	-r chr7:5566860 -m 10
	
echo "HANDLE NO READS IN INPUT"
java -Xmx500m -jar ../SamTextViewer.jar -i test_data/ds051.actb.bam -r chr7:5566860 -m 10 -f 16 -F 16

java -Xmx500m -jar ../SamTextViewer.jar -i test_data/ds051.actb.bam -r chr7:5566860 -m 10 -f 16

echo "CAN READ BAM WITH SORTED BY POSITION BUT MARKED AS SORTED BY queryname"
java -Xmx500m -jar ../SamTextViewer.jar -i test_data/ds051.unsorted_header.bam
