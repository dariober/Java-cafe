<!-- MarkdownTOC -->

- [SamTextViewer: A text based, no GUI genome viewer](#samtextviewer-a-text-based-no-gui-genome-viewer)
- [Usage](#usage)
  - [Quick start](#quick-start)
  - [Moving around the genome](#moving-around-the-genome)
  - [Display options](#display-options)
  - [Filtering reads](#filtering-reads)
  - [Searching features in annotation files](#searching-features-in-annotation-files)
  - [Note on regular expression](#note-on-regular-expression)
  - [Formatting of reads and features](#formatting-of-reads-and-features)
- [Getting help](#getting-help)
- [Supported input files](#supported-input-files)
- [Requirements and Installation](#requirements-and-installation)
  - [Installation quick start.](#installation-quick-start)
  - [A little more detail](#a-little-more-detail)
- [Performance](#performance)
- [TODO, FIXME etc](#todo-fixme-etc)
- [Credits](#credits)

<!-- /MarkdownTOC -->

# SamTextViewer: A text based, no GUI genome viewer


```SamTextViewer``` is a command line genome viewer useful to browse and visualize sequence alignment and annotation files
on **console screen**. It aims to be similar to ```samtools tview``` but with the flexibility of
GUI genome viewers like IGV.

Features that attempt to combine text based viewers (```tview```) with GUI viewers:

* Command line input and interaction, no graphical interface.
* Can load multiple files in various formats.
* Support for BS-Seq alignment files.
* Navigation and search options

![ex3](https://github.com/dariober/Java-cafe/blob/master/SamTextViewer/screenshots/ex3.png)

# Usage

## Quick start

These are some examples, for brevity using the helper `SamTextViewer`

Display a bam file together with a gtf annotation file, go straight to position chr7:5566640-5569055 (atcb gene). This is RNA-Seq data:

    SamTextViewer -r chr7:5566640-5569055 ds051.actb.bam hg19.gencode_genes_v19.gtf.gz

![ex1](https://github.com/dariober/Java-cafe/blob/master/SamTextViewer/screenshots/ex1.png)


The header line:

```
ds051.actb.bam; Each . = 44.68; max depth: 893.62x; 
```

gives the name of the file, scale of the read depth track (in this example one dot corresponds to 44.68 reads) and the maximum read depth
in the current view (893.62).

The line at the bottom of the tracks

```
chr7:5567688-5567847; 160 bp; 1.0 bp/char; Filters: -q 0 -f 0 -F 4; Mem: 552 MB;
```

shows the current position, the width and scale of the view, filters applied to the bam files and the memory usage.

For visualizing BS-Seq data add the `-bs` flag and provide a reference fasta file:

    SamTextViewer -bs fa chr7.fa ds051.actb.bam 


![exBSmode-2](https://github.com/dariober/Java-cafe/blob/master/SamTextViewer/screenshots/exBSmode-2.png)

![exBSmode](https://github.com/dariober/Java-cafe/blob/master/SamTextViewer/screenshots/exBSmode.png)

After starting `SamTextViewer` you can navigate the genome with the following interactive commands. Note that some options can be set either at start time or interactively, e.g. `-r`. 

## Moving around the genome

Use the following keys and commands to navigate the genome

```
f / b 
        Small step forward/backward 1/10 window
ff / bb
        Large step forward/backward 1/2 window
zi / zo
        Zoom in / zoom out
p / n
        Go to previous/next visited position
:<pos>
        Go to position <pos> on current chromosome
[+]/[-]<int>[k,m]
        Move forward/backward by <int> bases. Suffixes k and m allowed. E.g. -2m
```

To jump directly to a location use the -r option as `-r chr1:1000` or `chr1:1000-2000`

## Display options

These options can be used to set the limits of the y axis and their height.

```
visibile [show regex] [hide regex] [track regex]
        In annotation tracks, only include rows captured by [show regex] and exclude [hide regex].
        Apply to annotation tracks captured by [track regex]. With no optional arguments reset to default: "'.*' '^$' '.*'"
        Use '.*' to match everything and '^$' to hide nothing. Ex "visible .*exon.* .*CDS.* .*gtf#.*"
ylim <min> <max> [regex]
        Set limits of y axis for all track IDs captured by regex. Default regex: '.*'
dataCol <idx> [regex]
        Select data column for all bedgraph tracks captured by regex. <idx>: 1-based column index
```

also see:

```
print
        Turn on/off the printing of bed/gtf features in current interval
rNameOn / rNameOff
        Show/Hide read names
history
        Show visited positions
```

These options apply to bam files:

```
-m
    Maximum number of lines to print for read tracks. No limit If < 0
-rpm
    Toggle on/off the normalization of Reads Per Million for bam input. Default off
-d
    Maximum number of lines to print for coverage tracks
-ml
    Maximum number of lines to print for each methylation track
```

## Filtering reads

Reads in bam files can be filtered in and out using the -f, -F and -q options as in `samtools`:

```
-f
    Required sam flags. Use 4096 for reads on top strand
-F
    Filtering sam flags. Use 4096 for reads on top strand
-q
    Minumum mapping quality for a read to be considered

```

## Searching features in annotation files

These options apply to bed, gff and gtf files:

```
next <trackId>
        Move to the next feature in <trackId> on *current* chromosome
find <regex> [trackId]
        Find the next record in trackId matching regex. Use single quotes for strings containing spaces.
        For case insensitive matching prepend (?i) to regex e.g. '(?i).*actb.*'
```

Note that `find` searches the entire lines for matches to the given regular expression. So to find the 'ACTB' gene name use `find .*ACTB.*`, using just `find ACTB` as you would with `grep` will return no matches as no line in a bed file can match just that. 

## Note on regular expression

The options that take regular expressions assume some familiarity with regexes. In particular remember that lines are scanned as a single string for matches. These are some pointers:

* Almost always you want to surround your pattern with `.*`, e.g. to find the ACTB gene use `.*ACTB.*`, but note that this will hit LACTB as well unless you are more specific (e.g. `.*"ACTB".*`) or you use fancier regex. E.g. `.*[^A-Z0-9]ACTB[^A-Z0-9].*` will match ACTB but not LACTB or ACTB9. 

* Just using 'ACTB' as pattern, as you would do with `grep`, will return nothing as no line should contain only 'ACTB'.

* Use the `(?i)` modifier to match in case insensitve mode, e.g. '(?i).*actb.*'

* For reference, regexes are parsed by Java `String.match()` method.

## Formatting of reads and features

When aligned reads are show at single base resolution, read bases follow the same convention as samtools: 
Upper case letters and `.` for read align to forward strand, lower case and `,` otherwise; second-in-pair reads are underlined;
grey-shaded reads have mapping quality of <=5. In bisulfite mode the characters M, U, m, u are used for methylated and unmethylated bases on forward and reverse strands.

# Getting help

From command line:

```
SamTextViewer -h
```

Once `SamTextViewer` is started help can be displayed with `h`

```
SamTextViewer <input files>
...
chr1:1-160; 160 bp; 1.0 bp/char; Filters: -q 0 -f 0 -F 4; Mem: 633 MB; 
[h] for help: h   <-- h HERE to show help
```

# Supported input files

* **bam** and **cram** files should be sorted and indexed, e.g. with `samtools sort` and `samtools index`. Sam files are not supported
* **bedGraph** recognized by extension `.bedGraph` or `.bedgraph`
* **bigWig** recognized by extension `.bw` or `.bigWig`
* **bed**, **gtf**, **gff** recognized by respective extensions
* **tdf** This is very useful for quickly displaying very large intervals like tens of megabases or entire chromosomes see [tdf](https://www.broadinstitute.org/igv/TDF)
* Other extensions: Will be treated as bed files, provided the format is actually bed!

All plain text formats (bed, bedgraph, etc) can be read as gzipped and there is no need to decompress them.

Bedgraph files should be sorted by position a `sort -k1,1 -k2,2n` will do. Unindexed bedGraph files are first bgzipped and indexed to temporary files which are deleted on exit. This can take time 
for large files so consider creating the index once for all with [tabix](http://www.htslib.org/doc/tabix.html), *e.g.*

```
bgzip my.bedgraph &&
tabix -p bed my.bedgraph.gz
```

Bed & gtf file are not required to be sorted or index but in this case they are loaded in memory. To save memory and time for large files you can again index them as above. Loaded in memory is typically fast for files of up to ~1/2 million records.

For input format specs see also [UCSC format](https://genome.ucsc.edu/FAQ/FAQformat.html) and for guidelines on the choice of format see [IGV recommendations](https://www.broadinstitute.org/igv/RecommendedFileFormats).


# Requirements and Installation

## Installation quick start. 

In the commands below replace version number with the latest from [releases](https://github.com/dariober/Java-cafe/releases):

```
wget https://github.com/dariober/Java-cafe/releases/download/v0.1.0/SamTextViewer-0.1.0.zip
unzip SamTextViewer-0.1.0.zip
cp SamTextViewer.jar /usr/local/bin/ # Or ~/bin/
cp SamTextViewer /usr/local/bin/     # Or ~/bin/ 
```

## A little more detail

`SamTextViewer.jar` requires **Java 1.7+** and this is (should be) the only requirement. There is virtually no installation needed as `SamTextViewer` is pure Java and should work on most (all?) platforms. Download the zip file `SamTextViewer-x.x.x.zip` from [releases](https://github.com/dariober/Java-cafe/releases), unzip it and execute the jar file with

    java -jar /path/to/SamTextViewer.jar --help

To avoid typing `java -jar ...` every time, you can put both the helper 
script `SamTextViewer` and the jar file ```SamTextViewer.jar``` in the same directory in your `PATH` and execute with:

	SamTextViewer [options]

Note the helper is a bash script. To set the amount of memory available to java use the `-Xmx` option as e.g. `java -Xmx1500m -jar ...`.

If for some reason the text formatting misbehaves, disable it with the `-nf` option. 

# Performance

Alignment files are typically accessed very quickly but `SamTextViewer` becomes slow when the window size grows
above a few hundreds of kilobases. Annotation files (bed, gff, gtf) are loaded in memory unless they are indexed
with `tabix`. 

# TODO, FIXME etc

* Genomic coords are not checked against sam header!
* When setting `-d 0 -m 0` parsing bam files should be made faster  
* If a sam header is available:
  ** Add option to print it
  ** Add track showing the current position on the chromosome
* Header ` ylim: NaN, NaN; max: 884.69; .= 44.0;` can be confusing. E.g. if ylim is not set it show NaN.


# Credits

* Bam processing is mostly done with the [samtools/htsjdk](https://github.com/samtools/htsjdk) library.
* Bigwig and tdf are processed with classes from [IGV](https://github.com/igvteam/igv) source code.
* Block compression and indexing done using [jvarkit](https://github.com/lindenb/jvarkit)

