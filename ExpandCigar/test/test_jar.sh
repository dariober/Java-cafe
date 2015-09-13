#!/bin/bash

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

docstring="\nDESCRIPTION \n
Test the jar file works correctly \n
USAGE: \n
bash ./test_jar.sh /path/to/ExpandCigar.jar \n"

if [[ $1 == "" ]]
then
	echo -e $docstring
	exit 1
fi

path=`dirname $0`

# ------------------------------------------------------------------------------

echo ""
echo "CAN SHOW HELP"
java -jar -Xmx500m $1 -h

echo ""
echo "CAN SHOW VERSION"
java -jar -Xmx500m $1 -v

echo ""
echo "CAN PARSE BAM"
java -jar -Xmx500m $1 -i $path/test_data/aln.bam -o /dev/null
if [[ $? != 0 ]]
then
	echo "*** NON-ZERO STATUS ***"
	printf "${RED}FAILED ${NC}\n"
else 
	printf "${GREEN}PASSED ${NC}\n"
fi

echo ""
echo "CAN READ FROM STDIN"
samtools view -b $path/test_data/aln.bam | java -jar -Xmx500m $1 -i - -o /dev/null
if [[ $? != 0 ]]
then
	echo "*** NON-ZERO STATUS ***"
	printf "${RED}FAILED ${NC}\n"
else 
	printf "${GREEN}PASSED ${NC}\n"
fi

echo ""
echo "NM:i:0 GIVES NO CIGAR CONTAINING X"
out=`java -jar -Xmx500m $1 -i $path/test_data/aln.bam -o .sam | grep 'NM:i:0' | cut -f 6 | grep 'X'`
if [[ $out != "" ]]
then
	echo $out
	printf "${RED}FAILED ${NC}\n"
else 
	printf "${GREEN}PASSED ${NC}\n"
fi

echo ""
echo "SAME MPILEUP BEFORE AND AFTER EXPANSION"

# Remove every tag from bam. We want mpileup to work only with cigars and reads
samtools view -h test_data/aln.bam \
| cut -f1-11 \
| samtools view -Sb - \
| samtools mpileup -A - > /tmp/before.txt 

# Expand cigar
java -jar -Xmx500m $1 -i test_data/aln.bam | samtools mpileup -A - > /tmp/after.txt

echo "    * SAME POSITIONS AND DEPTH"
dfr=`diff <(cut -f1-4 /tmp/before.txt) <(cut -f1-4 /tmp/after.txt)`
if [[ $dfr != "" ]]
then
	echo "INCONSISTENCY FOUND"
	diff <(cut -f1-4 /tmp/before.txt) <(cut -f1-4 /tmp/after.txt)
	printf "${RED}FAILED ${NC}\n"
else 
	printf "${GREEN}PASSED ${NC}\n"
fi

echo "    * SAME CHAR COUNT"
cntBefore=`wc /tmp/before.txt | cut -d' ' -f1-3`
cntAfter=`wc /tmp/after.txt | cut -d' ' -f1-3`
if [[ $cntBefore != $cntAfter ]]
then
	echo "Different wc output"
	echo $cntBefore $cntAfter
	printf "${RED}FAILED ${NC}\n"
else 
	printf "${GREEN}PASSED ${NC}\n"
fi
rm /tmp/before.txt
rm /tmp/after.txt

echo ""
echo "CAN ADD PROG GROUP ENTRY IN HEADER"
pg=`java -jar -Xmx500m $1 -i test_data/aln.bam | samtools view -H - | grep '@PG' | grep 'ID:ExpandCigar'`

if [[ $pg == "" ]]
then
	echo "Program group not found"
	printf "${RED}FAILED ${NC}\n"
else 
	printf "${GREEN}PASSED ${NC}\n"
fi
