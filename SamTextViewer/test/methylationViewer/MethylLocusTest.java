package methylationViewer;

import static org.junit.Assert.*;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import samTextViewer.SamLocusIterator;
import samTextViewer.SamLocusIterator.LocusInfo;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Iterator;

import methylationViewer.MethylLocus;

import org.junit.Test;

import readWriteBAMUtils.ReadWriteBAMUtils;

public class MethylLocusTest {

	/*
	ds051.actb.bam; Max read depth: -1.0; Each . = -1.0x, -1.0 bp
	.............
	TT...........
	TT...........
	TT...........
	  ...........
	   ..........
	    .........
	    .........
	     ........
	           C.
	CATTTTTAAGGTG
	 5566781   5566791
	*/
	
	@Test
	public void canSetMethylLocus() {
		
		int from= 5566780;
		IndexedFastaSequenceFile faSeqFile = null;
		try {
			faSeqFile = new IndexedFastaSequenceFile(new File("test_data/chr7.fa"));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
				
		/* prepare locus info */
		SamReader samReader= ReadWriteBAMUtils.reader("test_data/ds051.actb.bam", ValidationStringency.SILENT);
		SAMFileHeader fh= samReader.getFileHeader();		
		IntervalList il= new IntervalList(fh);
		Interval interval= new Interval("chr7", from, from);
		il.add(interval);
		il.add(new Interval("chr7", 5566790, 5566790));
		SamLocusIterator samLocIter= new SamLocusIterator(samReader, il, true);
		Iterator<LocusInfo> iter= samLocIter.iterator();
		
		/* Test - See alignment above */ 
		LocusInfo locus= iter.next();
		MethylLocus m= new MethylLocus(locus, faSeqFile);
		assertEquals("C", m.getRefBase() + "");
		assertEquals(1, (int)m.getCntM());
		assertEquals(3, (int)m.getCntU());

		locus= iter.next();
		m= new MethylLocus(locus, faSeqFile);
		assertEquals("G", m.getRefBase() + "");
		assertEquals(0, (int)m.getCntM());
		assertEquals(0, (int)m.getCntU());		
	}
	
	@Test
	public void constructMethylLocusFromSubRefSeq(){
		
		SamReader samReader= ReadWriteBAMUtils.reader("test_data/ds051.actb.bam", ValidationStringency.SILENT);
		SAMFileHeader fh= samReader.getFileHeader();		
		IntervalList il= new IntervalList(fh);
		Interval interval= new Interval("chr7", 5566780, 5566780);
		il.add(interval);
		SamLocusIterator samLocIter= new SamLocusIterator(samReader, il, true);
		Iterator<LocusInfo> iter= samLocIter.iterator();
		
		LocusInfo locus = iter.next();
		byte[] subRefSeq= "CATTTTTAAGGTG".getBytes();
		int offset= 5566780;
		MethylLocus m= new MethylLocus(locus, subRefSeq, offset);	
		assertEquals("C", m.getRefBase() + "");
		assertEquals(1, (int)m.getCntM());
		assertEquals(3, (int)m.getCntU());
	}
	
}
