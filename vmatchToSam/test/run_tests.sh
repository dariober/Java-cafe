#!/bin/bash

docstring="\nTest runner for vmatchToSam.jar. Requires picard in ~/bin/picard.jar
\nUsage:
\ncd /.../Java-cafe/trunk/vmatchToSam
\nbash test/run_test.sh /path/to/vmatchToSam.jar \n"

if [[ $1 == "" ]];
then
  echo -e $docstring
  exit
fi

echo -e "\nCAN SHOW HELP"
java -jar $1 -h

echo -e "\nCAN CONVERT VMATCH TO SAM"
java -jar $1 -fa test_data/ref.fa test_data/aln_1.vmatch.txt > test.tmp.sam
java -jar ~/bin/picard.jar ValidateSamFile IGNORE=MISSING_PLATFORM_VALUE I=test.tmp.sam 2> /dev/null

echo -e "\nCAN CONVERT FROM STDIN"
cat  test_data/aln_1.vmatch.txt | java -jar $1 -fa test_data/ref.fa > test.tmp.sam
java -jar ~/bin/picard.jar ValidateSamFile IGNORE=MISSING_PLATFORM_VALUE I=test.tmp.sam 2> /dev/null

rm test.tmp.sam