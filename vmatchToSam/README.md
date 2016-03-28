## Utility to convert vmatch output to sam format

[vmatch](http://www.vmatch.de/) is a versatile software tool for eﬃciently solving large scale sequence matching tasks.
`vmatchToSam` converts the vmatch output to [sam format](http://samtools.github.io/hts-specs/SAMv1.pdf) so that
various tools for sam analysis can be used on vmatch.

## Input and Output

Currently, the input for `vmatchToSam` is the vmatch output in sequence format (`-s`) and with the `-showdesc` option enabled.
`vmatchToSam` accepts stdin and writes to stdout.

For example a typical use case might be:

```
vmatch -s -showdesc 0 -p -d -complete -q query.fa ref.fa \
| java -jar ~/bin/vmatchToSam.jar -fa ref.fa > out.vmatch.sam
```

#### vmatch options explained:

* `-s` print the sequence alignment. This option is currently required.

* `-showdesc 0` show reference sequence names. Default is to show the index position of the sequence in the reference file.
This option is not strictly necessary but if not present `vmatchToSam.jar` will not be able to add the correct header, unless an *ad hoc*
list of reference names and lengths is used. I don't see much point in not using the `-showdesc` option.

* `-p -d` align query sequences as they are (`-d`, direct) and as reverse complement (`-p`, palindromic).
Not required, but you almost always want to enable both these options. 

* `-complete` Force the query sequences to align completely to the reference.
At least one of the options `-complete`, `-hxdrop`, `-exdrop`, or `-l` is required by `vmatch`.

* `-q reads.fa` Query sequences in fasta format.

* `ref.fa` Reference sequence indexed using the vmatch tool `mkvtree`, typical executed as
`mkvtree -dna -allout -pl -db ref.fa`.

`vmatch` options not compatible with `vmatchToSam` are: `-nodist`, `-noevalue`, `-noscore`, `-noidentity` and formats other then the sequence alignment.

#### vmatchToSam options and output

* `-fa ref.fa` Reference sequences, used to generate the sam header. Alternatively the `-fai ref.fa.fai` option can be used where
`ref.fa.fai` is a tab separated file of sequence names and lengths. Options `-fa` and `-fai` are not required but if not given the sam
output will be headerless.


The sam output will have mapping quality 255 (unavailable); the cigar string uses `=` and `X` for match/mismatch instead of the generic `M`.
The following attributes are set:

* `XE:f` E-value

* `NM:i` Edit distance

* `XP:f` Percent match

* `AS:i` Alignment score

* `RG:Z:NA` Read group always set to NA.

See also [addAlignmentTagsToBam.py](https://github.com/dariober/bioinformatics-cafe/blob/master/addAlignmentTagsToBam.py) for a tool to add some useful tags to sam.

## Install and Requirements

`vmatchToSam` is written in Java so it should work out of the box. Download the compiled jar file from
https://github.com/dariober/Java-cafe/releases and execute as:

```
wget https://github.com/dariober/Java-cafe/releases/download/vx.x.x/vmatchToSam.jar
java -jar /path/to/vmatchToSam.jar -h
```

## TODO

* Rank alignments within queries by alignment score or else so it's easy to get the "best" query aligmnent.

* Accept xml output from vmatch. Problem is xml does not contain the aligned sequences.