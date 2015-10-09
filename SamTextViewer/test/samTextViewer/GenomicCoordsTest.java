package samTextViewer;

import static org.junit.Assert.*;

import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Locale;
import java.util.Random;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.ValidationStringency;

import org.junit.Test;

import readWriteBAMUtils.ReadWriteBAMUtils;

public class GenomicCoordsTest {
	
	public static SAMSequenceDictionary samSeqDict= ReadWriteBAMUtils
			.reader("test_data/ds051.short.bam", ValidationStringency.STRICT)
			.getFileHeader().getSequenceDictionary();
	
	@Test
	public void canZoom(){
		GenomicCoords gc= new GenomicCoords("chr1:101-105", samSeqDict); 
		gc.zoomOut();
		assertEquals(95, (int)gc.getFrom());
		assertEquals(111, (int)gc.getTo());
		for(int i= 0; i < 30; i++){
			gc.zoomOut();
		}
		assertEquals(1, (int)gc.getFrom());
		assertEquals(samSeqDict.getSequence("chr1").getSequenceLength(), (int)gc.getTo()); // Doesn't extend beyond chrom
		
		gc= new GenomicCoords("chr1:100-200", samSeqDict);
		gc.zoomIn();
		assertEquals(125, (int)gc.getFrom());
		assertEquals(175, (int)gc.getTo());
		
		for(int i= 0; i < 20; i++){
			gc.zoomIn();
		}
		assertEquals(149, (int)gc.getFrom());
		assertEquals(151, (int)gc.getTo());
		

		
	}
	
	@Test
	public void canGoToRegionAndReset(){
		GenomicCoords gc= GenomicCoords.goToRegion("chr1:1-10", "test_data/ds051.short.bam", 100);
		assertEquals("chr1", gc.getChrom());
		assertEquals(1, (int)gc.getFrom());
		assertEquals(10, (int)gc.getTo());
		
		gc= GenomicCoords.goToRegion("chr1:1", "test_data/ds051.short.bam", 100);
		assertEquals("chr1", gc.getChrom());
		assertEquals(1, (int)gc.getFrom());
		assertEquals(100, (int)gc.getTo());
		
		gc= GenomicCoords.goToRegion("chr1:1", "test_data/ds051.short.bam", 100);
		assertEquals("chr1", gc.getChrom());
		assertEquals(1, (int)gc.getFrom());
		assertEquals(100, (int)gc.getTo());

		gc= GenomicCoords.goToRegion("chr1:100000-100200", "test_data/ds051.short.bam", 100);
		
		assertEquals("chr1", gc.getChrom());
		assertEquals(100000, (int)gc.getFrom());
		assertEquals(100200, (int)gc.getTo());

		gc= GenomicCoords.goToRegion("chrM:16550", "test_data/ds051.short.bam", 100);
		assertEquals(16472, (int)gc.getFrom());
		assertEquals(16571, (int)gc.getTo());

		gc= GenomicCoords.goToRegion("chrM:16550-17000", "test_data/ds051.short.bam", 100);
		assertEquals(16472, (int)gc.getFrom());
		assertEquals(16571, (int)gc.getTo());

		gc= GenomicCoords.goToRegion("nonsense", "test_data/ds051.short.bam", 100);
		assertEquals(null, gc.getChrom());
		assertEquals(null, gc.getFrom());
		assertEquals(null, gc.getTo());

//		gc= GenomicCoords.goToRegion("chrM:-1000", "test/test_data/ds051.short.bam", 100);
//		System.out.println(gc);
//		assertEquals(null, gc.getChrom());
//		assertEquals(null, gc.getFrom());
//		assertEquals(null, gc.getTo());
		
	}
	
	@Test 
	public void returnNullForChromNotFoundInSam(){
		
		GenomicCoords gc= new GenomicCoords("nonsense:1-10", samSeqDict);
		assertEquals(null, gc.getChrom());
		assertEquals(null, gc.getFrom());
		assertEquals(null, gc.getTo());
		
		gc= new GenomicCoords("nonsense", samSeqDict);
		assertEquals(null, gc.getChrom());
		assertEquals(null, gc.getFrom());
		assertEquals(null, gc.getTo());	
	}

	@Test
	public void canParseAndGetCoords() {		
		
		String region= "chrX:1-100";
		GenomicCoords gc= new GenomicCoords(region, samSeqDict);
		assertEquals("chrX", gc.getChrom());
		assertEquals(1, (int)gc.getFrom());
		assertEquals(100, (int)gc.getTo());
		
		region= "chrX : 1 -    100";
		
		gc= new GenomicCoords(region, samSeqDict);
		assertEquals("chrX", gc.getChrom());
		assertEquals(1, (int)gc.getFrom());
		assertEquals(100, (int)gc.getTo());
		
		region= "chrX:1,000,000-1,100,000";
		gc= new GenomicCoords(region, samSeqDict);
		assertEquals(1000000, (int)gc.getFrom());
		assertEquals(1100000, (int)gc.getTo());

		region= "chrX:1,000,000-1,100,000";
		gc= new GenomicCoords(region, samSeqDict);
		assertEquals(1000000, (int)gc.getFrom());
		
		region= "chrX:-1";
		gc= new GenomicCoords(region, samSeqDict);
		//assertEquals(null, gc.getFrom());

	}
	
	@Test
	public void canHandleNonSenseCoords(){
		
		// These cases should be user input errors.
		
		String region= "chrX:1,0,0,0,000-1100000"; // Take care, this should fail
		GenomicCoords gc= new GenomicCoords(region, samSeqDict);
		assertEquals(1000000, (int)gc.getFrom());

		region= "chrX:110-100"; // start > end
		gc= new GenomicCoords(region, samSeqDict);
		assertEquals(null, gc.getFrom());

		region= "chrX:-1-1000"; // Negative start: Non sense
		gc= new GenomicCoords(region, samSeqDict);
		assertEquals(null, gc.getFrom());

		region= "chrX:0-1000"; // 
		gc= new GenomicCoords(region, samSeqDict);
		assertEquals(null, gc.getFrom());

		region= "chrX: 1,100,000,000"; // 
		gc= new GenomicCoords(region, samSeqDict);
		assertEquals(155270560, (int)gc.getFrom()); // Reset to chrom size
		assertEquals(null, gc.getTo());
	}
	
	@Test
	public void canHandleChromsWithMetachars(){
	
		String weird= "chr:X-Y";
		samSeqDict.addSequence(new SAMSequenceRecord(weird, 10000));
		String region= weird + ":1-100"; // Names with colon and hyphen ok.
		GenomicCoords gc= new GenomicCoords(region, samSeqDict);
		assertEquals(weird, gc.getChrom());
		assertEquals(1, (int)gc.getFrom());
		assertEquals(100, (int)gc.getTo());
	}
	
	@Test
	public void canGetChromAndFromWithoutTo(){

		String region= " chrX : 100 ";
		GenomicCoords gc= new GenomicCoords(region, samSeqDict);
		assertEquals("chrX", gc.getChrom());
		assertEquals(100, (int)gc.getFrom());
		assertEquals(null, gc.getTo());
		
	}
	
	@Test
	public void canParseStringToGenomicCoords(){
		String x= "\"bla:192121-10\":10-100";
		x= "chr7:10-100";
		x= "chr7:1-1000";
		GenomicCoords gc= new GenomicCoords(x, samSeqDict);	
	}
	
	@Test
	public void canGetChromOnly(){
		String region= " chrX";
		GenomicCoords gc= new GenomicCoords(region, samSeqDict);
		assertEquals("chrX", gc.getChrom());
		assertEquals(null, gc.getFrom());
		assertEquals(null, gc.getTo());		
	}


	@Test//(expected=NumberFormatException.class)
	public void formatNotValid(){ // This should be better handled. What if chrom names contain ':'?
		String region= "chrX:chr2";
		GenomicCoords gc= new GenomicCoords(region, samSeqDict);
		assertEquals("chrX", gc.getChrom());
		assertEquals(null, gc.getFrom());
		assertEquals(null, gc.getTo());
	}
	
	@Test
	public void canResetCoords(){
		String region= "chrM:1000000";
		GenomicCoords gc= new GenomicCoords(region, samSeqDict);
		gc.correctCoordsAgainstSeqDict(samSeqDict);
		assertEquals(null, gc.getTo());
		assertEquals(16571, (int)gc.getFrom());
		
		region= "chrM:20000-30000";
		gc= new GenomicCoords(region, samSeqDict);
		gc.correctCoordsAgainstSeqDict(samSeqDict);
		assertEquals(16571, (int)gc.getTo());
		assertEquals(16571, (int)gc.getFrom());

		region= "nonsense";
		gc= new GenomicCoords(region, samSeqDict);
		gc.correctCoordsAgainstSeqDict(samSeqDict);
		assertEquals(null, gc.getChrom());
		assertEquals(null, gc.getTo());
		assertEquals(null, gc.getFrom());

	}
		
}
