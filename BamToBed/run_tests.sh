#!/bin/bash

# diff commands should return nothing.

echo "CAN CONVERT ALL RECORDS EXCEPT UNMAPPED"

n_exp=`samtools view -b test_data/grm056_test.bam | bamToBed | wc -l`
n_obs=`java -jar ~/Tritume/BamToBed.jar -i test_data/grm056_test.bam | wc -l`
diff <(echo $n_exp) <(echo $n_obs)


echo "CAN PARSE BY CHROM AND FLAGS"

samtools view -b -q 41 -F 16 test_data/grm056_test.bam 10:20099950 | bamToBed > expected.tmp.bed
java -jar ~/Tritume/BamToBed.jar -q 41 -F 16 -chrom 10 -from 20099950 -i test_data/grm056_test.bam > observed.tmp.bed
diff  expected.tmp.bed observed.tmp.bed
rm expected.tmp.bed observed.tmp.bed