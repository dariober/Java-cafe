package samTextViewer;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import tracks.Track;

import org.junit.Test;

import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;

public class UtilsTest {

	static SamReaderFactory srf=SamReaderFactory.make();
	static SamReader samReader= srf.open(new File("test_data/ds051.short.bam"));
	public static SAMSequenceDictionary samSeqDict= samReader.getFileHeader().getSequenceDictionary();
	
	public static String fastaFile= "test_data/chr7.fa";

	@Test
	public void canSetYlimits() throws InvalidCommandLineException{
		String cmdInput= ".*.gz$ 0 10 :ylim";
		List<Track> tracks= new ArrayList<Track>();
		Track t1= new Track(); t1.setFilename("foo.gz"); tracks.add(t1);
		Track t2= new Track(); t2.setFilename("foo.txt"); tracks.add(t2);
		Track t3= new Track(); t3.setFilename("bla.gz"); tracks.add(t3);

		Utils.setTrackYlimitsForRegex(cmdInput, tracks);
		
		assertEquals(0, tracks.get(0).getYmin(), 0.001);
		assertEquals(10, tracks.get(0).getYmax(), 0.001);
		assertEquals(Double.NaN, tracks.get(1).getYmax(), 0.001);
	}
	
	@Test
	public void canGetFileTypeFromName(){
		assertEquals("bigWig",
		Utils.getFileTypeFromName("/Users/berald01/Downloads/wgEncodeCaltechRnaSeqGm12878R2x75Il400SigRep2V2.bigWig"));
	} 

	@Test
	public void canInitRegion() throws IOException{
		assertEquals("chrM", Utils.initRegionFromFile("test_data/ds051.short.bam"));
		assertEquals("chr1", Utils.initRegionFromFile("test_data/refSeq.hg19.short.bed"));
		assertEquals("chr1", Utils.initRegionFromFile("test_data/refSeq.hg19.short.sort.bed.gz"));
		assertEquals("chr9", Utils.initRegionFromFile("test_data/hg18_var_sample.wig.v2.1.30.tdf"));
		assertEquals("chr1", Utils.initRegionFromFile("/Users/berald01/Downloads/wgEncodeCaltechRnaSeqGm12878R2x75Il400SigRep2V2.bigWig"));		
	}
	
	@Test
	public void testBamHasIndex() throws IOException{
		assertTrue(Utils.bamHasIndex("test_data/ds051.short.bam"));
		assertTrue(!Utils.bamHasIndex("test_data/ds051.noindex.bam"));
	}
	
	@Test
	public void canGetClosestIndex(){
		int windowSize= 150;
		List<Double> mapping = Utils.seqFromToLenOut(1, 1000000, windowSize);
		for(int i=0; i < windowSize; i++){
			System.out.println("Index: " + i + " position: " + mapping.get(i));
		}
		
		assertEquals(windowSize-1, Utils.getIndexOfclosestValue(1000000, mapping));
		assertEquals((windowSize/2)-1, Utils.getIndexOfclosestValue(1000000/2, mapping));
		assertEquals(0, Utils.getIndexOfclosestValue(1, mapping));
		assertEquals(0, Utils.getIndexOfclosestValue((6712/2.0)-1, mapping));
		assertEquals(1, Utils.getIndexOfclosestValue((6712/2.0)+1, mapping));
		assertEquals(0, Utils.getIndexOfclosestValue((6712/2.0), mapping));
		assertEquals(102, Utils.getIndexOfclosestValue(684564, mapping));
		assertEquals(112, Utils.getIndexOfclosestValue(750000, mapping));
	}
	
	@Test
	public void canGenerateSequence(){
		Utils.seqFromToLenOut(15, 17, 13);
		Utils.seqFromToLenOut(17, 15, 13);
		Utils.seqFromToLenOut(15, 26, 13).size();
		Utils.seqFromToLenOut(15, 15, 13);
		
		// Length of 1 returns "from" like R seq(0, 10, length.out= 1) -> 0
		assertEquals(1, Utils.seqFromToLenOut(0, 10, 1).size());
		assertEquals(0, Utils.seqFromToLenOut(0, 10, 1).get(0), 0.00001);		
	}
	
	@Test
	public void canParseInputAndUpdateGenomicCoords() throws InvalidGenomicCoordsException, IOException{
		GenomicCoords gc= new GenomicCoords("chr7:100-200", samSeqDict, 100, fastaFile);

		//String region= Utils.parseConsoleInput("-r chr8:1-1000", gc);
		//assertEquals("chr8:1-1000", region);
		
		//String region= Utils.parseConsoleInput("-r chr8:1", gc);
		//assertEquals("chr8:1", region);

		//region= Utils.parseConsoleInput("-r chr8", gc);
		//assertEquals("chr8", region);
		
		String region= Utils.parseConsoleInput("10", gc);
		assertEquals("chr7:110-210", region);

		region= Utils.parseConsoleInput("-1000", gc);
		// assertEquals("chr7:110-210", region);
		
		String rawInput= "-r chr10 -F 1024";
		List<String> clArgs= Arrays.asList(rawInput.split("\\s+"));
		// System.out.println(clArgs.indexOf("-R"));		
	}
	
	// @Test
	//public void canParseColonOperator(){
	//	// Test : operator
	//	GenomicCoords gc= new GenomicCoords("chr7:100-200", samSeqDict);
	//	assertEquals("chr7:200", Utils.parseConsoleInput(":200", gc));
	//}
	//
	// @Test
	//public void canparseConsoleInputOperator_f_b(){
		// Test : operator
	//	GenomicCoords gc= new GenomicCoords("chr7:100-200", samSeqDict);
	//	assertEquals("chr7:110-210", Utils.parseConsoleInput("f", gc));
	//	assertEquals("chr7:90-190", Utils.parseConsoleInput("b", gc));
	//	// Smallest step is 1bp
	//	gc= new GenomicCoords("chr7:100-102", samSeqDict);
	//	assertEquals("chr7:101-103", Utils.parseConsoleInput("f", gc));
	//	assertEquals("chr7:99-101", Utils.parseConsoleInput("b", gc));
	//}
	
	//public void canGetStartOfBamFile(){
	//	GenomicCoords gc= Utils.getStartCoordsOfBAM("test_data/ds051.short.bam");
	//	assertEquals("chr7", gc.getChrom());
	//	assertEquals(5566778, (int)gc.getFrom());
	//	assertEquals(gc.getFrom(), gc.getTo());
	//}
	//
	//public void canGetStartOfChrom(){
	//	GenomicCoords gc= Utils.getStartCoordsOfBAM("test_data/ds051.short.bam", "chr7");
	//	assertEquals("chr7", gc.getChrom());
	//	assertEquals(5566778, (int)gc.getFrom());
	//	assertEquals(gc.getFrom(), gc.getTo());
	//	
	//	// null for chrom with no reads. Fails if chrom does not exist in header.
	//	gc= Utils.getStartCoordsOfBAM("test_data/ds051.short.bam", "chrY");
	//	assertEquals(null, gc.getChrom());
	//	assertEquals(null, gc.getFrom());
	//	assertEquals(null, gc.getTo());
	//}
	

}