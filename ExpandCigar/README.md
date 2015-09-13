# Description

This program expands the 'M' operator in cigar string to the operators **X** (aligned with mismatch) and **=** (aligned with match). 
To work correctly it is necessary that the MD tag and read sequences are present and correct. 
If the MD tag is absent it can be calculated with `samtools calmd` and the output piped to expand the cigar. For example:

```
samtools calmd aln.bam ref.fa | java -Xmx500m -jar ExpandCigar.jar -i - > out.bam
```

The runnable jar can be found in [releases](https://github.com/dariober/Java-cafe/releases).
