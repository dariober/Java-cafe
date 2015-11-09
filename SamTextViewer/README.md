# SamTextViewer 


```SamTextViewer``` is a command line genome viewer to visualize sequence alignment and annotation files
on console screen. It aims to be similar to ```samtools tview``` but with the flexibility of
GUI genome viewers like IGV.

Features that attempt to combine text based viewers (```tview```) with GUI viewers:

* Command line input and interaction, no graphical interface.
* Can load multiple files in various formats.
* Support for BS-Seq alignent files.
* Navigation and search options

![ex3](https://github.com/dariober/Java-cafe/blob/master/SamTextViewer/screenshots/ex3.png)

# Installation

There is virtually no installation needed as ```SamTextViewer.jar``` is pure Java. 
Download the zip file `SamTextViewer-x.x.x.zip` from [releases](https://github.com/dariober/Java-cafe/releases), unzip it and execute the jar file with

    java -jar /path/to/SamTextViewer.jar --help

To avoid typing ```java -jar ...``` every time, you can put both the helper 
file `SamTextViewer` and the jar file ```SamTextViewer.jar``` in a directory in your PATH and execute with:

	SamTextViewer [options]

# Usage 

These are some examples, for brevity using the helper `SamTextViewer`

View help

    SamTextViewer -h

Display a bam file together with a gtf annotation file, go straight to position chr7:5566640-5569055 (atcb gene). This is RNA-Seq data:

    SamTextViewer -r chr7:5566640-5569055 ds051.actb.bam hg19.gencode_genes_v19.gtf.gz

![ex1](https://github.com/dariober/Java-cafe/blob/master/SamTextViewer/screenshots/ex1.png)


The header line:

```
ds051.actb.bam; Each . = 44.68; max depth: 893.62x; 
```

gives the name of the file, scale of the read depth track (in this example one dot coresponds to 44.68 reads) and the maximum read depth
in the current view (893.62).

The line at the bottom of the tracks

```
chr7:5567688-5567847; 160 bp; 1.0 bp/char; Filters: -q 0 -f 0 -F 4; Mem: 552 MB;
```

shows the current position, the width and scale of the view, filters applied to the bam files and the memory usage.

For visualizing BS-Seq data add the `-bs` flag and provide a reference fasta file:

    SamTextViewer -bs fa chr7.fa ds051.actb.bam 

After starting `SamTextViewer` you can navigate the genome with the following interactive commands. Note that some options can be set either at start time
or interctively, e.g. `-r`.

    Navigation options

    [f]/[b]
            Small, step forward/backward 1/10 window
    [ff]/[bb]
            Large step forward/backward 1/2 window
    [zi]/[zo]
            Zoom in / zoom out
    [p]/[n]
           Go to previous/next visited position
    [:pos]
            Go to position <pos> on current chromosome
    +/-<int>
            Move forward/backward by <int> bases. Suffixes k and m allowed. E.g. -2m
    <filename> [:n]
            Move to the next feature in <filename> on current chromosome
    <filename> <string> [:find]
            Find the next record in file containing string. Use single quotes for strings containing spaces.
    [print]
            Turn on/off the printing of bed/gtf features in current interval
    [rNameOn]/[rNameOff]
            Show/Hide read names
    [history]
            Show visited positions
    [q]
            Quit
    [h]
            Show this help
    
    Command line options
    
    -r
        Go to region. Format 1-based as 'chrom:start-end' or 'chrom:start' or 'chrom'
    -f
        Required sam flags. Use 4096 for reads on top strand
    -F
        Filtering sam flags. Use 4096 for reads on top strand
    -q
        Minumum mapping quality for a read to be considered
    -m
        Maximum number of lines to print for read tracks. No limit If < 0
    -rpm
        Toggle on/off the normalization of Reads Per Million for bam input. Default off
    -d
        Maximum number of lines to print for coverage tracks. No limit if < 0
    -ml
        Maximum number of lines to print for each methylation track. No limit if < 0


# Supported input files

For input format specs see also [UCSC format](https://genome.ucsc.edu/FAQ/FAQformat.html) and for choice of format see [IGV recommendations](https://www.broadinstitute.org/igv/RecommendedFileFormats)

* **bam** and **cram** files should be sorted and indexed, e.g. with `samtools sort` and `samtools index`. Sam files are not supported.
* **bedGraph** recognized by extension `.bedGraph`, can be uncompressed or gzipped.
* **bigWig** recognized by extension `.bw` or `.bigWig`
* **bed**, **gtf**, **gff** recognized by respective extensions, can be gzipped.
* **tdf** This is very useful for quickly displaying very large intervals (tens of megabases).
* Other extensions: Will be treated as bed files, provided the format is actually bed!

Large bed, gtf, bedGraph files should be sorted, bgzipped and indexed with **`tabix`** for fast access and memory efficiency
(see [tabix manual](http://www.htslib.org/doc/tabix.html)). If not indexed, these files will be loaded in memory, which is anyway fine
for files of up to ~1/2 million records.

# Performance

Alignment files are typically accessed very quickly but `SamTextViewer` becomes slow when the displayed window size grows
above a few hundreds of kilobases. Consider setting `-d 0` and `-m 0` to temporarily turn off the visualization of bam files.

# Credits

Bam processing is mostly done with the [samtools/htsjdk](https://github.com/samtools/htsjdk) library.
Bigwig and tdf are processed with classes from [IGV](https://github.com/igvteam/igv) source code.

