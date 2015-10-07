# SamTextViewer 

```SamTextViewer``` is a command line genome viewer designed to visualize seuence alignment 
on console screen. It aims to be similar to ```samtools tview``` but with the flexibility of
GUI genome viewers like IGV. ```SamTextViewer``` allows to quickly visualize bam files on 
remote server without the need of GUI.

Features that attempt to combine text based viewers (```tview```) with GUI viewers:

* Command line input and interaction. Fully text based. 
* Can load multiple bam files to be visualized in parallel 
* Support for BS-Seq alignent files.

# Installation

There is virtually no installation needed as ```SamTextViewer.jar``` is written in Java. 
Download the compile jar file from [releases](https://github.com/dariober/Java-cafe/releases) and execute with

    java -jar /path/to/SamTextViewer.jar --help

To avoid typing ```java -jar ...``` every time, you can put the helper 
file [SamTextViewer](https://github.com/dariober/Java-cafe/SamTextViewer/) and the jar
 file ```SamTextViewer.jar``` in direcotry in your PATH and execute with:

	SamTextViewer --help

# Usage 

Input bam files must sorted and indexed, e.g. with ```samtools```. The fasta file of the reference genome 
must be as well indexed, e.g. with ```samtools faidx myref.fa```.

## Loading bam and reference

View multiple bam files in parallel

    SamTextViewer -i aln1.bam aln2.bam [alnN.bam]

By default the viewer will go to the start of the first bam file, like ```samtools view aln.bam```. 
To go directly to a given position use the ```-r``` option as ```chrom``` or ```chrom:pos```

   SamTextViewer -i aln1.bam -r chr18:10000

## Read filtering

```SamTextViewer```

 