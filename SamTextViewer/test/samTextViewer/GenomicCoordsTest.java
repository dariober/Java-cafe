package samTextViewer;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.junit.Test;

import exceptions.InvalidGenomicCoordsException;

public class GenomicCoordsTest {
	
	static SamReaderFactory srf=SamReaderFactory.make();
	static SamReader samReader= srf.open(new File("test_data/ds051.short.bam"));
	public static SAMSequenceDictionary samSeqDict= samReader.getFileHeader().getSequenceDictionary();
	
	public static String fastaFile= "test_data/chr7.fa";
	
	@Test
	public void canPrintChromMap() throws InvalidGenomicCoordsException, IOException{
			
		GenomicCoords gc= new GenomicCoords("chr7", 1, 1, samSeqDict, 10, null);
		
		String chromMap= gc.getChromIdeogram();		
		assertEquals("*--------|", chromMap);
		
		
		gc= new GenomicCoords("chr7", 1, 1000000000, samSeqDict, 10, null);
		chromMap= gc.getChromIdeogram();		
		assertEquals("**********", chromMap);
		
		
		gc= new GenomicCoords("chr7", 200000000, 200000000, samSeqDict, 10, null);
		chromMap= gc.getChromIdeogram();		
		assertEquals("|--------*", chromMap);
		
		
		gc= new GenomicCoords("chr7", 20000000, 55000000, samSeqDict, 16, null);
		chromMap= gc.getChromIdeogram();		
		assertEquals("|-****---------|", chromMap);
		
	}
	
	@Test
	public void printRefSeq() throws InvalidGenomicCoordsException, IOException{
		GenomicCoords gc= new GenomicCoords("chr7", 5540580, 5540590, null, 100, "test_data/chr7.fa");
		//GenomicCoords gc= new GenomicCoords("chr7", 5540580, 5540590, samSeqDict, 100, null);
		System.out.println("START");
		System.out.println(gc.printableRefSeq(true));
		System.out.println("DONE");
	}
	
	@Test
	public void canTestForEqualCoords() throws InvalidGenomicCoordsException, IOException{
		GenomicCoords gc= new GenomicCoords("chr1", 1, 10, null, 100, null);
		GenomicCoords other= new GenomicCoords("chr1", 1, 10, null, 1000, null);
		assertTrue(gc.equalCoords(other));
		
		GenomicCoords other2= new GenomicCoords("chr2", 1, 10, null, 1000, null);
		assertTrue(!gc.equalCoords(other2));
		other2= new GenomicCoords("chr1", 2, 10, null, 1000, null);
		assertTrue(!gc.equalCoords(other2));
		other2= new GenomicCoords("chr1", 1, 100, null, 1000, null);
		assertTrue(!gc.equalCoords(other2));
		
		GenomicCoords gc2= (GenomicCoords) gc.clone();
		System.out.println(gc);
		System.out.println(gc2);
		gc.zoomOut();
		System.out.println(gc);
		System.out.println(gc2);
		// assertTrue(gc2.equalCoords(gc));		
	}
	
	@Test
	public void canInitializeSamSeqDictFromGenomeFile() throws IOException{
	
		List<String> insam= new ArrayList<String>();
		// From resource:
		assertEquals(93, GenomicCoords.getSamSeqDictFromAnyFile(insam, null, "hg19").size());
		// From bam header:
		assertEquals(25, GenomicCoords.getSamSeqDictFromAnyFile(insam, null, "test_data/ds051.short.bam").size());
	}
	
	@Test
	public void canGetSamSeqDict() throws IOException{
		List<String> insam= new ArrayList<String>();
		insam.add("test_data/ds051.short.bam.bai"); // This will not produce anything
		insam.add("test_data/ds051.short.bam");
		SAMSequenceDictionary ssd = GenomicCoords.getSamSeqDictFromAnyFile(insam, null, null);
		assertEquals(25, ssd.size());
		
		// From indexed fasta
		insam= new ArrayList<String>();
		insam.add("test_data/ds051.short.bam.bai"); // This will not produce anything
		ssd = GenomicCoords.getSamSeqDictFromAnyFile(null, fastaFile, null);
		assertEquals(1, ssd.size());		
	}
	
	@Test
	public void canPrintRefSeq() throws InvalidGenomicCoordsException, IOException{
		GenomicCoords gc= new GenomicCoords("chr7", 5566770, 5566790, samSeqDict, 100, fastaFile);
		assertEquals("CACTTGGCCTCATTTTTAAGG\n", gc.printableRefSeq(true));
		// with format
		gc= new GenomicCoords("chr7", 5566770, 5566772, samSeqDict, 100, fastaFile);
		assertEquals("[107;31mC[0m[107;34mA[0m[107;31mC[0m\n", gc.printableRefSeq(false));
	}
	
	@Test
	public void canConstructGenomicCoords() throws InvalidGenomicCoordsException, IOException{
		
		GenomicCoords gc= new GenomicCoords("chr7", 1, 100, samSeqDict, 1000, fastaFile);
		assertEquals("chr7", gc.getChrom());
		assertEquals(1, (int)gc.getFrom());
		assertEquals(100, (int)gc.getTo());

		gc= new GenomicCoords("chr7:1-100", samSeqDict, 1000, fastaFile);
		assertEquals("chr7", gc.getChrom());
		assertEquals(1, (int)gc.getFrom());
		assertEquals(100, (int)gc.getTo());

		gc= new GenomicCoords("chr7", 11, null, samSeqDict, 1000, fastaFile);		
		assertEquals(1010, (int)gc.getTo());

		gc= new GenomicCoords("chr7", null, null, samSeqDict, 1000, fastaFile);		
		assertEquals(1, (int)gc.getFrom());
		assertEquals(1000, (int)gc.getTo());

		gc= new GenomicCoords("chr7", 1000000000, 1000000000, samSeqDict, 1000, fastaFile); // Reset to size of chrom
		assertEquals(159138663, (int)gc.getFrom());

		gc= new GenomicCoords("chr7:100", samSeqDict, 1000, fastaFile);
		assertEquals(100, (int)gc.getFrom());
		assertEquals(100+1000-1, (int)gc.getTo());

		gc= new GenomicCoords("chr7:1,000,000,000", samSeqDict, 1000, null);
		assertEquals(159138663-1000+1, (int)gc.getFrom()); // Reset to chrom size
		assertEquals(159138663, (int)gc.getTo());

		gc= new GenomicCoords("chr7:1,000,000,000", null, 1000, null);
		assertEquals(1000000000, (int)gc.getFrom()); // Fine, no dict to check against.
		assertEquals(1000000000+1000-1, (int)gc.getTo());

		
	}
	
	@Test
	public void canGetRefSeq() throws InvalidGenomicCoordsException, IOException{
		GenomicCoords gc= new GenomicCoords("chr7", 5566770, 5566790, samSeqDict, 1000, fastaFile);
		assertEquals("CACTTGGCCTCATTTTTAAGG", new String(gc.getRefSeq()));
		gc= new GenomicCoords("chr7", 5566770, 5566790, samSeqDict, 20, fastaFile);
		// System.out.println(gc.getBpPerScreenColumn());
		assertEquals(null, gc.getRefSeq());
	}
	
	@Test(expected = InvalidGenomicCoordsException.class)
	public void canThrowNullChrom() throws InvalidGenomicCoordsException, IOException {
		new GenomicCoords(null, 1, 100, samSeqDict, 100, fastaFile);
	}
	
	@Test(expected = InvalidGenomicCoordsException.class)
	public void canThrowInvalidCoords() throws InvalidGenomicCoordsException, IOException {
		new GenomicCoords("chr7", -1, 100, samSeqDict, 100, fastaFile);
	}
	
	@Test(expected = InvalidGenomicCoordsException.class)
	public void canThrowChromNotInDict() throws InvalidGenomicCoordsException, IOException {
		new GenomicCoords("nonsense", 1, 100, samSeqDict, 100, fastaFile);
	}
	
	@Test
	public void canZoom() throws InvalidGenomicCoordsException, IOException{
		GenomicCoords gc= new GenomicCoords("chr1:101-105", samSeqDict, 100, null); 
		gc.zoomOut();
		assertEquals(99, (int)gc.getFrom()); // exp 95,111 if zoom fact is x2
		assertEquals(107, (int)gc.getTo()); // 
		for(int i= 0; i < 40; i++){
			gc.zoomOut();
		}
		assertEquals(1, (int)gc.getFrom());
		assertEquals(samSeqDict.getSequence("chr1").getSequenceLength(), (int)gc.getTo()); // Doesn't extend beyond chrom
		
		gc= new GenomicCoords("chr1:101-1000", samSeqDict, 100, null);
		gc.zoomIn();
		assertEquals(326, (int)gc.getFrom());
		assertEquals(776, (int)gc.getTo());
		
		// Zoom-in in small interval has no effect
		gc= new GenomicCoords("chr1:101-200", samSeqDict, 200, null);
		gc.zoomIn();
		assertEquals(101, (int)gc.getFrom());
		assertEquals(200, (int)gc.getTo());
		
		gc= new GenomicCoords("chr1:1-200", samSeqDict, 200, null);
		gc.zoomIn();
		assertEquals(1, (int)gc.getFrom());
		assertEquals(200, (int)gc.getTo());

		gc= new GenomicCoords("chrM:16561-16571", samSeqDict, 200, null); // End of chrom
		gc.zoomIn();
		assertEquals(16561, (int)gc.getFrom());
		assertEquals(16571, (int)gc.getTo());
	}
	
	@Test
	public void canPrepareRuler() throws InvalidGenomicCoordsException, IOException{
		GenomicCoords gc= new GenomicCoords("chr1:101-110", samSeqDict, 7, null);
		assertEquals(7, gc.getMapping().size());
		assertEquals(101.0, gc.getMapping().get(0), 0.01);
		assertEquals(102.5, gc.getMapping().get(1), 0.01);
		assertEquals(102.5, gc.getMapping().get(1), 0.01);
		assertEquals(1.42857142, gc.getBpPerScreenColumn(), 0.01);
	}
	
	@Test
	public void canPrintRuler() throws InvalidGenomicCoordsException, IOException{
		
		GenomicCoords gc= new GenomicCoords("chr1:101-200", samSeqDict, 50, null);
		assertEquals(50, gc.printableRuler(10).length());
	
	}
	
	/*
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
		
		GenomicCoords gc= new GenomicCoords("nonsense:1-10", samSeqDict, 100);
		assertEquals(null, gc.getChrom());
		assertEquals(null, gc.getFrom());
		assertEquals(null, gc.getTo());
		
		gc= new GenomicCoords("nonsense", samSeqDict, 100);
		assertEquals(null, gc.getChrom());
		assertEquals(null, gc.getFrom());
		assertEquals(null, gc.getTo());	
	}

	@Test
	public void canParseAndGetCoords() {		
		
		String region= "chrX:1-100";
		GenomicCoords gc= new GenomicCoords(region, samSeqDict, 100);
		assertEquals("chrX", gc.getChrom());
		assertEquals(1, (int)gc.getFrom());
		assertEquals(100, (int)gc.getTo());
		
		region= "chrX : 1 -    100";
		
		gc= new GenomicCoords(region, samSeqDict, 100);
		assertEquals("chrX", gc.getChrom());
		assertEquals(1, (int)gc.getFrom());
		assertEquals(100, (int)gc.getTo());
		
		region= "chrX:1,000,000-1,100,000";
		gc= new GenomicCoords(region, samSeqDict, 100);
		assertEquals(1000000, (int)gc.getFrom());
		assertEquals(1100000, (int)gc.getTo());

		region= "chrX:1,000,000-1,100,000";
		gc= new GenomicCoords(region, samSeqDict, 100);
		assertEquals(1000000, (int)gc.getFrom());
		
		region= "chrX:-1";
		gc= new GenomicCoords(region, samSeqDict, 100);
		//assertEquals(null, gc.getFrom());

	}
	
	@Test
	public void canHandleNonSenseCoords(){
		
		// These cases should be user input errors.
		
		String region= "chrX:1,0,0,0,000-1100000"; // Take care, this should fail
		GenomicCoords gc= new GenomicCoords(region, samSeqDict, 100);
		assertEquals(1000000, (int)gc.getFrom());

		region= "chrX:110-100"; // start > end
		//gc= new GenomicCoords(region, samSeqDict, 100);
		//assertEquals(null, gc.getFrom());

		region= "chrX:-1-1000"; // Negative start: Non sense
		gc= new GenomicCoords(region, samSeqDict, 100);
		assertEquals(null, gc.getFrom());

		region= "chrX:0-1000"; // 
		gc= new GenomicCoords(region, samSeqDict, 100);
		assertEquals(null, gc.getFrom());

		region= "chrX: 1,100,000,000"; // 
		gc= new GenomicCoords(region, samSeqDict, 100);
		assertEquals(155270560, (int)gc.getFrom()); // Reset to chrom size
		assertEquals(null, gc.getTo());
	}
	
	@Test
	public void canHandleChromsWithMetachars(){
	
		String weird= "chr:X-Y";
		samSeqDict.addSequence(new SAMSequenceRecord(weird, 10000));
		String region= weird + ":1-100"; // Names with colon and hyphen ok.
		GenomicCoords gc= new GenomicCoords(region, samSeqDict, 100);
		assertEquals(weird, gc.getChrom());
		assertEquals(1, (int)gc.getFrom());
		assertEquals(100, (int)gc.getTo());
	}
	
	@Test
	public void canGetChromAndFromWithoutTo(){

		String region= " chrX : 100 ";
		GenomicCoords gc= new GenomicCoords(region, samSeqDict, 100);
		assertEquals("chrX", gc.getChrom());
		assertEquals(100, (int)gc.getFrom());
		assertEquals(null, gc.getTo());
		
	}
	
	@Test
	public void canParseStringToGenomicCoords(){
		String x= "\"bla:192121-10\":10-100";
		x= "chr7:10-100";
		x= "chr7:1-1000";
		GenomicCoords gc= new GenomicCoords(x, samSeqDict, 100);	
	}
	
	@Test
	public void canGetChromOnly(){
		String region= " chrX";
		GenomicCoords gc= new GenomicCoords(region, samSeqDict, 100);
		assertEquals("chrX", gc.getChrom());
		assertEquals(null, gc.getFrom());
		assertEquals(null, gc.getTo());		
	}


	@Test//(expected=NumberFormatException.class)
	public void formatNotValid(){ // This should be better handled. What if chrom names contain ':'?
		String region= "chrX:chr2";
		GenomicCoords gc= new GenomicCoords(region, samSeqDict, 100);
		assertEquals("chrX", gc.getChrom());
		assertEquals(null, gc.getFrom());
		assertEquals(null, gc.getTo());
	}
	
	@Test
	public void canResetCoords(){
		String region= "chrM:1000000";
		GenomicCoords gc= new GenomicCoords(region, samSeqDict, 100);
		gc.correctCoordsAgainstSeqDict(samSeqDict);
		assertEquals(null, gc.getTo());
		assertEquals(16571, (int)gc.getFrom());
		
		region= "chrM:20000-30000";
		gc= new GenomicCoords(region, samSeqDict, 100);
		gc.correctCoordsAgainstSeqDict(samSeqDict);
		assertEquals(16571, (int)gc.getTo());
		assertEquals(16571, (int)gc.getFrom());

		region= "nonsense";
		gc= new GenomicCoords(region, samSeqDict, 100);
		gc.correctCoordsAgainstSeqDict(samSeqDict);
		assertEquals(null, gc.getChrom());
		assertEquals(null, gc.getTo());
		assertEquals(null, gc.getFrom());

	} */
			
}
