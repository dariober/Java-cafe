# Mark duplicates by read 5' and 3'end position

`MarkDupsByStartEnd.jar` marks *single end* reads as duplicates if they share the
same **start and end** position. 
In contrast, `Picard/MrkDuplicates` and `samtools rmdup` consider SE reads as duplicates if their 5'end
are the same and they ignore the 3'end. This approach is overly conservative if read length is comparable
to or larger then the library fragment size.

**NB** Paired reads are returned to output unchanged. 
Also reads marked as unmapped, failed, supplementary are not considered.

## Recommended usage

Reads should be trimmed to remove adapter sequences before alignment. 
However, low quality 3'ends should not be trimmed as they would artificially make the sequence fragment size
shorter. `bwa mem` copes well with low quality ends.

By default reads are duplicates if they share: chrom, start, end, strand, read group. Read group can be ignored with `-rg`.

**NB** If the same library is sequenced more than once with different read length it is not correct to consider the 3'end. 
Use `Picard/MarkDuplicates` in this case. 

Usage:

    java -jar MarkDupsByStartEnd.jar -h
    java -jar MarkDupsByStartEnd.jar -i <aln.bam|aln.sam> -o <md.bam|md.sam>
    samtools view -u aln.bam | java -jar MarkDupsByStartEnd.jar -i - | samtools view ...

You probably want to use `java -Xmx2g`. 
Input need not to be sorted. Output will be coordinate sorted by default.

## Requirements 

GNU Coreutils `sort` is required to be on `PATH`. `sort` is used to sort reads after unclipping.
Use `sort` version 8.6+ to take advantage of parallel threading.

# TODO

* `sort`: Use `getOutputStream` to feed sort instead of using a tmp file?
* ~~Picard MarkDuplicates use the sum of base qualities to pick best alignment. 
This is probably the right thing to do since reads with genuine SNPs would be penalized.~~
* ~~Read from stdin~~. Memo: Use `samtools view -u aln.bam | ...` (or `-b`) instead of simply `samtools view aln.bam | ...` 
otherwise the output bam will be headerless. 
* ~~Function to get tmp file name in the same dir as the input file. Like `<input>.xxxx.tmp`. *I.e.* make sure the tmp file does not exists already!~~
* ~~Clean up logging etc.~~
 