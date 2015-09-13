# Description

This program expands the 'M' operator in cigar string to the operators **X** (aligned with mismatch) and **=** (aligned with match). 
To work correctly it is necessary that the MD tag and read sequences are present and correct. 
If the MD tag is absent it can be calculated with `samtools calmd` and the output piped to expand the cigar. For example:

```
samtools calmd aln.bam ref.fa | java -Xmx500m -jar ExpandCigar.jar -i - > out.bam
```

The runnable jar can be found in [releases](https://github.com/dariober/Java-cafe/releases).

Program's help (might be outdated):
```
usage: ExpandCigar [-h] --insam INSAM [--outsam OUTSAM]
                   [--oldCigarTag OLDCIGARTAG] [--version]

DESCRIPTION
Expand cigar string in input bam file to convert M operator to X and =.

For source code and further information see:
https://github.com/dariober/Java-cafe/tree/master/ExpandCigar


optional arguments:
  -h, --help             show this help message and exit
  --insam INSAM, -i INSAM
                         Input sam or bam file. Use - to read from stdin.
  --outsam OUTSAM, -o OUTSAM
                         Output file. Format will  be  sam or bam depending
                         on extension.
                         Use - to print BAM to  stdout or '.sam' for sam to
                         stdout. (default: -)
  --oldCigarTag OLDCIGARTAG, -t OLDCIGARTAG
                         Store the original cigar  in  this tag, e.g. 'XX'.
                         If this tag exists it will be overwritten.
  --version, -v

```
