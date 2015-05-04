# Mark duplicates by read 5' and 3'end position

Get compiled `MarkDupsByStartEnd.jar` from [releases](https://github.com/dariober/Java-cafe/releases).

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
    
    samtools view -u aln.bam \
    | java -jar MarkDupsByStartEnd.jar -i - \
    | samtools view ...

If reading from stdin remember to include the sam header, e.g. `samtools view -u/-b/-h`, otherwise
the output will be headerless. 

You probably want to use `java -Xmx2g`. 
Input need not to be sorted. Output will be coordinate sorted by default.

## Help

(Might be outdated).

	usage: MarkDupsByStartEnd [-h] --insam INSAM [--outsam OUTSAM] [--unsortedOutput] 
							  [--ignoreReadGroup] [--validationStringency {SILENT,LENIENT,STRICT}]
	                          [--version]
	
	DESCRIPTION
	Mark duplicate reads by looking at both the unclipped read start and end position.
	For further details see
	https://github.com/dariober/Java-cafe/tree/master/MarkDupsByStartEnd
	
	optional arguments:
	  -h, --help             show this help message and exit
	  --insam INSAM, -i INSAM
	                         Input sam or bam file. Use - to read from stdin.
	  --outsam OUTSAM, -o OUTSAM
	                         Output file.
	                         Format will be sam or bam depending on extension.
	                         Use - to print SAM to stdout. (default: -)
	  --unsortedOutput, -us  If set, output reads will be unsorted.
	                         By default reads are coordinate sorted. (default: false)
	  --ignoreReadGroup, -rg
	                         Ignore read group info. If set, positional duplicates sharing *different* 
	                         read groups will be considered duplicates. (default: false)
	  --validationStringency {SILENT,LENIENT,STRICT}, -vs {SILENT,LENIENT,STRICT}
	                         Set picard validation stringency level. (default: SILENT)
	  --version, -v


## TODO

* Add arg to set max number of records in memory/ Profile memory usage vs performance for increasing memory.
