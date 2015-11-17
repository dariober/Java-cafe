package tracks;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.junit.Test;

import exceptions.InvalidGenomicCoordsException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import samTextViewer.GenomicCoords;
import samTextViewer.SamLocusIterator;
import tracks.TrackCoverage;

public class TrackCoverageTest {
	
	public static List<SamRecordFilter> filters= new ArrayList<SamRecordFilter>();
	
	static SamReaderFactory srf=SamReaderFactory.make();
	static SamReader sr= srf.open(new File("test_data/ds051.short.bam"));
	public static SAMSequenceDictionary samSeqDict= sr.getFileHeader().getSequenceDictionary();
	public static String fastaFile= "test_data/chr7.fa";
	
	//@Test
	public void testSpeedSamLocIter() throws IOException{
		
		SamReader samReader= srf.open(new File("/Volumes/My_Passport_for_Mac/tmp/rhh147-148_untreat_14102014_atac_hacat.bam"));
		SAMFileHeader fh= samReader.getFileHeader();

		// FAST
		// IntervalList il= new IntervalList(fh);
		// int x= 5000000;
		// while(x < 10000000){
		//	il.add(new Interval("chr7", x, x));
		//	x += 10;
		// }
		//IntervalList il= new IntervalList(fh);
		// SamLocusIterator samLocIter= new SamLocusIterator(samReader, il); // SLOW
		//System.out.println(il.size());

		IntervalList il= new IntervalList(fh);
		il.add(new Interval("chr7", 5000000, 60000000));
		SamLocusIterator samLocIter= new SamLocusIterator(samReader, il); // FAST
		System.out.println("SamLocIter DONE");
		
		long t0= System.currentTimeMillis();
		Iterator<samTextViewer.SamLocusIterator.LocusInfo> iter= samLocIter.iterator(); // FAST
		long t1= System.currentTimeMillis();
		System.out.println("Iterator done in: " + (t1- t0));
		
		// FAST
		int i= 0;
		long nbp= 0;
		while(iter.hasNext()){
			i++;
			samTextViewer.SamLocusIterator.LocusInfo locusInfo= iter.next();
			nbp += locusInfo.getRecordAndPositions().size();
		}
		long t2= System.currentTimeMillis();
		System.out.println(i + " loci Done in: " + (t2- t1) + " ms; counted " + nbp);
		samLocIter.close();
		samReader.close();
	}
	
	@Test
	public void testRPMnorm() throws InvalidGenomicCoordsException, IOException{
		
		int yMaxLines= 11;
		int windowSize= 101;

		GenomicCoords gc= new GenomicCoords("chr7", 5566770, 5566870, samSeqDict, windowSize, fastaFile);
		TrackCoverage tc= new TrackCoverage("test_data/ds051.short.bam", gc, filters, false);
		tc.setyMaxLines(yMaxLines);
		tc.setRpm(true);
		assertEquals(1000000, tc.getMaxDepth(), 0.1);
	}
	
	@Test
	public void canPrintCoverageTrack() throws IOException, InvalidGenomicCoordsException {
		int yMaxLines= 11;
		int windowSize= 101;

		GenomicCoords gc= new GenomicCoords("chr7", 5566770, 5566870, samSeqDict, windowSize, fastaFile);
		TrackCoverage tc= new TrackCoverage("test_data/ds051.short.bam", gc, filters, false);
		tc.setyMaxLines(yMaxLines);
		System.out.println(gc.toString());
		System.out.println(tc.printToScreen());
		System.out.println(tc.getScorePerDot());
	}

	@Test 
	public void testWithZeroReadsRegion() throws InvalidGenomicCoordsException, IOException{
		
		GenomicCoords  gc= new GenomicCoords("chr7", 1, 1000, samSeqDict, 20, null);
		TrackCoverage tc= new TrackCoverage("test_data/ds051.short.bam", gc, filters, false);
		tc.setyMaxLines(10);
		assertEquals("____________________", tc.printToScreen());
	}

	@Test
	public void canPrintCoverage() throws IOException, InvalidGenomicCoordsException {
		int yMaxLines= 11;
		int windowSize= 101;

		GenomicCoords gc= new GenomicCoords("chr7:5568363-5568390", samSeqDict, windowSize, fastaFile);
		TrackCoverage tc= new TrackCoverage("test_data/ds051.actb.bam", gc, filters, false);
		tc.setyMaxLines(yMaxLines);
		System.out.println(gc.toString());
		System.out.println(tc.printToScreen());
		System.out.println(tc.getScorePerDot());
	}
}
