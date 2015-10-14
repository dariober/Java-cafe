package coverageViewer;

import static org.junit.Assert.*;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import org.apache.commons.lang3.StringUtils;
import org.junit.Test;

import readWriteBAMUtils.ReadWriteBAMUtils;
import samTextViewer.GenomicCoords;
import samTextViewer.SamLocusIterator;
import samTextViewer.SamLocusIterator.LocusInfo;

public class CoverageViewerTest {

	public static List<SamRecordFilter> filters= new ArrayList<SamRecordFilter>();
	
	public static SAMSequenceDictionary samSeqDict= ReadWriteBAMUtils
			.reader("test_data/ds051.short.bam", ValidationStringency.STRICT)
			.getFileHeader().getSequenceDictionary();
	
	
	@Test
	public void stubCovergeViewer(){
		
		String sam= "test_data/ds051.short.bam";
		String chrom= "chr7";
		int from= 5566781;
		int to= from+10;
		
		GenomicCoords gc= new GenomicCoords(chrom, from, to, samSeqDict);
		
		CoverageViewer cw= new CoverageViewer(sam, gc, 1000, filters);
		assertEquals(to - from + 1, cw.getDepth().size());
		assertEquals("[4.0, 5.0, 6.0, 8.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 10.0]", cw.getDepth().toString());
		assertEquals(10, cw.getMaxDepth(), 0.001);
		assertEquals(5566781, (int)cw.getGenomicPositions().get(0));
		
		String expProfile= 
"    ......:\n" +
"   ::::::::\n" +
" .:::::::::\n" +
":::::::::::\n" +
":::::::::::";

		String obsProfile= StringUtils.join(cw.getProfileStrings(10), "\n");
		assertEquals(expProfile, obsProfile);

		expProfile= ""+
"   .:::::::\n" +
":::::::::::";		
		obsProfile= StringUtils.join(cw.getProfileStrings(2), "\n");
		assertEquals(expProfile, obsProfile);
		
	}
	
	@Test
	public void canCompressCoverageView(){
		
		String sam= "test_data/ds051.short.bam";
		String chrom= "chr7";
		int from= 5566781;
		int to= from+200;
		int windowSize= 10;
		GenomicCoords gc= new GenomicCoords(chrom, from, to, samSeqDict);
		
		CoverageViewer cw= new CoverageViewer(sam, gc, windowSize, filters);
		cw.getMappingToScreen(gc.getMapping(windowSize));
		cw.getProfileStrings(10);
	}
	
	//@Test
	public void canConstructLargeIntervalBySampling(){
		
		String sam= "test_data/ds051.actb.bam";
		String chrom= "chr7";
		int from= 5566778;
		int to= 5567778;
		int windowSize= 50;
		GenomicCoords gc= new GenomicCoords(chrom, from, to, samSeqDict);
		CoverageViewer cw= new CoverageViewer(sam, gc, windowSize, filters);
		// The number of loci collected is ~ LOC_PER_WINDOW * windowSize
		assertEquals(CoverageViewer.LOC_PER_WINDOW * windowSize, cw.getDepth().size(), 200);
		
		// Start sampling at approx start of coords
		assertEquals(from, cw.getGenomicPositions().get(0), 50);
		// Stop sampling at approx end of coords.
		assertEquals(to, cw.getGenomicPositions().get(cw.getGenomicPositions().size()-1), 50);
	}
	
	@Test
	public void canCompressLargeCoverageView(){
		
		String sam= "test_data/ds051.actb.bam";
		String chrom= "chr7";
		int from= 5566781;
		int to= from+10000;
		int windowSize= 50;
		GenomicCoords gc= new GenomicCoords(chrom, from, to, samSeqDict);
		
		CoverageViewer cw= new CoverageViewer(sam, gc, windowSize, filters);
		cw.getMappingToScreen(gc.getMapping(windowSize));
		assertEquals(from, cw.getGenomicPositions().get(0), 500);
		assertEquals(to, cw.getGenomicPositions().get(cw.getGenomicPositions().size()-1), 500);
		
		System.out.println(cw.getMappingToScreen(gc.getMapping(windowSize)));
		System.out.println(cw.getDepth().size());

		/* ================ 
		SamReader samReader= ReadWriteBAMUtils.reader("test_data/long.sort.bam", ValidationStringency.SILENT);
		SAMFileHeader fh= samReader.getFileHeader();
		IntervalList il= new IntervalList(fh);
		il.add(new Interval("chr7", 1, 5602271));
		SamLocusIterator samLocIter= new SamLocusIterator(samReader, il, true);
		long t0= System.currentTimeMillis();
		int n= 0;
		while(samLocIter.hasNext()){
			LocusInfo x = samLocIter.next();
			x.getRecordAndPositions().size();
			n++;
		}
		long t1= System.currentTimeMillis();
		System.out.println(t1-t0 + " " + n);
		*/
	}	
}
