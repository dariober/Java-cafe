package methylationViewer;

import static org.junit.Assert.*;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import samTextViewer.SamLocusIterator;
import samTextViewer.SamLocusIterator.LocusInfo;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import methylationViewer.MethylLoci;
import methylationViewer.MethylLocus;

import org.apache.commons.lang3.StringUtils;
import org.junit.Test;

import coverageViewer.CoverageViewer;
import readWriteBAMUtils.ReadWriteBAMUtils;
import samTextViewer.Utils;

public class MethylLociTest {

	public static List<SamRecordFilter> filters= new ArrayList<SamRecordFilter>();
	
	@Test
	public void canConstructFromCoverageViewer() {
		String sam= "test_data/ds051.actb.bam";
		String chrom= "chr7";
		int from= 5566780;
		int to= from+9;
		int windowSize= 100;
		
		CoverageViewer cw= new CoverageViewer(sam, chrom, from, to, windowSize, filters);
		byte[] subRefSeq= "CATTTTTAAGGTG".getBytes();
		int offset= 5566780;
		MethylLoci ml= new MethylLoci(cw, subRefSeq, offset);
		//System.out.println(ml);
		assertEquals(1.0, (float)ml.getMDepth().get(0), 0.1);
		assertEquals(3.0, (float)ml.getUDepth().get(0), 0.1);
		assertEquals(to-from+1, ml.getRefBases().size());
	}
	
	@Test
	public void canConstructFromMethylLociList() {

		byte[] subRefSeq= "CATTTTTAAgGTG".getBytes();
	
		SamReader samReader= ReadWriteBAMUtils.reader("test_data/ds051.actb.bam", ValidationStringency.SILENT);
		SAMFileHeader fh= samReader.getFileHeader();		
		IntervalList il= new IntervalList(fh);
		il.add(new Interval("chr7", 5566780, 5566780));
		il.add(new Interval("chr7", 5566781, 5566781));
		il.add(new Interval("chr7", 5566782, 5566782));
		il.add(new Interval("chr7", 5566783, 5566783));
		il.add(new Interval("chr7", 5566784, 5566784));
		il.add(new Interval("chr7", 5566785, 5566785));
		il.add(new Interval("chr7", 5566786, 5566786));
		il.add(new Interval("chr7", 5566787, 5566787));
		il.add(new Interval("chr7", 5566788, 5566788));
		il.add(new Interval("chr7", 5566789, 5566789));
		il.add(new Interval("chr7", 5566790, 5566790));
		il.add(new Interval("chr7", 5566791, 5566791));
		il.add(new Interval("chr7", 5566792, 5566792));
		SamLocusIterator samLocIter= new SamLocusIterator(samReader, il, true);
		Iterator<LocusInfo> iter= samLocIter.iterator();
		
		List<MethylLocus> mlist= new ArrayList<MethylLocus>();
		while(iter.hasNext()){
			MethylLocus m= new MethylLocus(iter.next(), subRefSeq, 5566780);
			mlist.add(m);
		}
		
		/*Construct */
		MethylLoci ml= new MethylLoci(mlist);
		assertEquals(5566780, (int)ml.getGenomicPositions().get(0));
		//System.out.println(ml);
	}

	@Test
	public void canPrintMethylationProfileString() {

		byte[] subRefSeq= "CATTTTTAAGGTG".getBytes();

		String sam= "test_data/ds051.actb.bam";
		String chrom= "chr7";
		int from= 5566780;
		int to= from+subRefSeq.length - 1;
		int windowSize= 5;
		
		CoverageViewer cw= new CoverageViewer(sam, chrom, from, to, windowSize, filters);
		int offset= 5566780;
		MethylLoci ml= new MethylLoci(cw, subRefSeq, offset);
		//System.out.println(ml);
		String obsProfile= StringUtils.join(ml.getMethylProfileStrings(10, true), "\n");
		//System.out.println(obsProfile);
		
		//ArrayList<Float> xl= new ArrayList<Float>();
		//xl.add(null); xl.add(null); xl.add((float)2); xl.add((float)4); 
		//Float x= Utils.calculateAverage(xl);
		// System.out.println( x );
	}

	@Test
	public void canCompressMethylLoci() {

		byte[] subRefSeq= "CATTTTTAAGGTG".getBytes();

		String sam= "test_data/ds051.actb.bam";
		String chrom= "chr7";
		int from= 5566780;
		int to= from+subRefSeq.length - 1;
		int windowSize= 5;
		
		CoverageViewer cw= new CoverageViewer(sam, chrom, from, to, windowSize, filters);
		int offset= 5566780;
		MethylLoci ml= new MethylLoci(cw, subRefSeq, offset);
		System.out.println(ml);
		ml.compressCovergeViewer(windowSize);
		System.out.println(ml);
		String obsProfile= StringUtils.join(ml.getMethylProfileStrings(10, true), "\n");
		//System.out.println(obsProfile);
	
		int x= (int)Math.round(Float.NaN);
		System.out.println(x == 0);
		
		//ArrayList<Float> xl= new ArrayList<Float>();
		//xl.add(null); xl.add(null); xl.add((float)2); xl.add((float)4); 
		//Float x= Utils.calculateAverage(xl);
		// System.out.println( x );
	}

	
}
