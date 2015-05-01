#!/bin/bash

set -e

pathToJar=$1
if [[ $1 == "" ]]; 
then
echo -e "\nThis is a test suite"
echo -e "\nUSAGE"
echo "markDupsbyStartEndTester.sh /path/to/MarkDupsByStartEnd.jar"
echo ""
exit 1
fi

shopt -s expand_aliases
alias markdups='java -jar '$1' -vs STRICT'

## START TESTING
## =============

echo "
0: CAN SHOW VERSION"
markdups --version

## ALL diff should return empty

echo "
1: NOTHING MARKED OR CHANGED => INPUT==OUTPUT"
# --------------------------------------------------
markdups -i testfile-1.sam > out.sam
out=`diff testfile-1.sam out.sam`
if [[ $out == '' ]]; then 
	echo "PASSED"
else
	echo "******** FAILED **********"; exit
fi
rm out.sam

echo "
2: IGNORING READ GROUP"
# ---------------------------
markdups -i testfile-1.sam -rg > out.sam
out=`diff expected-1.sam out.sam`
if [[ $out == '' ]]; then 
	echo "PASSED"
else
	echo "******** FAILED **********"; exit
fi
rm out.sam

echo "
3: OUTPUT IS INDEXABLE HENCE SORTED"
# ----------------------------------
markdups -i testfile-1.sam -rg -o out.bam && 
samtools index out.bam &&
rm out.bam out.bam.bai

echo "
4: CAN READ FROM STDIN"
# ---------------------
samtools view -bS testfile-1.sam | markdups -i - | samtools view -S - > /dev/null
 
echo "
5: SOFT CLIPPING IS EXTENDED"
# ---------------------------
markdups -i testfile-2.sam > out.sam
out=`diff expected-2.sam out.sam `
if [[ $out == '' ]]; then 
	echo "PASSED"
else
	echo "******** FAILED **********"; exit
fi
rm out.sam

echo "
6: READS ARE SCORED BY BASE QUALITY THEN BY MAPQ"
# -----------------------------------------------
markdups -i testfile-3.sam > out.sam
out=`diff expected-3.sam out.sam`
if [[ $out == '' ]]; then 
	echo "PASSED"
else
	echo "******** FAILED **********"; exit
fi
rm out.sam

echo "
7: PE, FAILED, UNMAPPED, SUPP READS ARE IGNORED"
# ----------------------------------------------
markdups -i testfile-4.sam > out.sam
out=`diff testfile-4.sam out.sam`
if [[ $out == '' ]]; then 
	echo "PASSED"
else
	echo "******** FAILED **********"; exit
fi
rm out.sam

echo "
8: SE READS ARE UN-MARKED"
# ------------------------
markdups -i testfile-5.sam > out.sam
out=`diff expected-5.sam out.sam`
if [[ $out == '' ]]; then 
	echo "PASSED"
else
	echo "******** FAILED **********"; exit
fi
rm out.sam
