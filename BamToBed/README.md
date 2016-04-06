## Utility to convert BAM file to BED

Convert input BAM file to BED in a similar way as [bedtools bamtobed](http://bedtools.readthedocs.org/en/latest/content/tools/bamtobed.html). 

Features not available in `bedtools bamtobed`: 

* Allows random access of chromosomes and positions (similar to `samtools view ... [reg]`) 
* Allows filtering for mapq, and bitwise flags (similar to `samtools view -q <int> -f <int> -F <int>`)
* No installation required as this is a Java program

Output is bed format with columns: chrom, start, end, read name, mapq, strand.

Example usage: Keep reads with mapq >= 5, discard second in pair reads, extract chr18 only:

```
java -jar BamToBed.jar -i in.bam -q 5 -F 128 -chrom chr18 \
| <pipe to other prog>
```

See also [explain sam flag](https://broadinstitute.github.io/picard/explain-flags.html)

## Help 

_Might be outdated_

```
java -jar  BamToBed.jar -h
usage: BamToBed [-h] [--inbam INBAM] [--chrom CHROM] [--from FROM] [--to TO] [--mapq MAPQ] [--requiredFlag REQUIREDFLAG] [--filterFlag FILTERFLAG] [--version]

DESCRIPTION
Convert BAM file to BED
See also https://github.com/dariober/Java-cafe/tree/master/BamToBed


optional arguments:
  -h, --help             show this help message and exit
  --inbam INBAM, -i INBAM
                         Input bam file. Should be sorted and indexed (default: -)
  --chrom CHROM, -chrom CHROM
                         Only output reads on this chromosome (default: )
  --from FROM, -from FROM
                         Only output reads starting from this position. 1-based, 0 for start of reference (default: 0)
  --to TO, -to TO        Only output reads up to this position. 1-based, 0 for end of reference (default: 0)
  --mapq MAPQ, -q MAPQ   Minimum mapq score to output a read. Equivalent to samtools view -q (default: 0)
  --requiredFlag REQUIREDFLAG, -f REQUIREDFLAG
                         Output reads with these bits set in flag. Equivalent to `samtools view -f` (default: 0)
  --filterFlag FILTERFLAG, -F FILTERFLAG
                         Do not output reads with these bits set in flag. Equivalent to `samtools view -F`. Unmapped reads are always discared. (default: 0)
  --version, -v
```

## Installation and Requirement

No installation required as this is pure java. JVM 1.6+ required. Just download the compiled [BamToBed.jar](https://github.com/dariober/Java-cafe/releases/download/v0.1.0/BamToBed.jar) from Releases and execute as

```
java -jar ~/path/to/BamToBed.jar ...
```
