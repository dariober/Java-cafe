Start project

# TODO

* `sort`: Use `getOutputStream` to feed sort instead of using a tmp file?
* ~~Picard MarkDuplicates use the sum of base qualities to pick best alignment. 
This is probably the right thing to do since reads with genuine SNPs would be penalized.~~
* ~~Read from stdin~~. Memo: Use `samtools view -u aln.bam | ...` (or `-b`) instead of simply `samtools view aln.bam | ...` 
otherwise the output bam will be headerless. 
* ~~Function to get tmp file name in the same dir as the input file. Like `<input>.xxxx.tmp`. *I.e.* make sure the tmp file does not exists already!~~
* ~~Clean up logging etc.~~
 