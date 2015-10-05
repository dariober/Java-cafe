package samTextViewer;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Set;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.SamLocusIterator;
import htsjdk.samtools.util.SamLocusIterator.LocusInfo;
import htsjdk.samtools.util.SamLocusIterator.RecordAndOffset;

import org.junit.Test;

import readWriteBAMUtils.ReadWriteBAMUtils;

public class UtilsTest {

	public static SAMSequenceDictionary samSeqDict= ReadWriteBAMUtils
			.reader("test/test_data/ds051.short.bam", ValidationStringency.STRICT)
			.getFileHeader().getSequenceDictionary();
	
	@Test
	public void testStackReads() {
		
		/*
NNNNNNNNNNNNNNNNNNNN
AAAAA
 CCCCC
      TTTTT
        GGGGG
              AAAAA
		 */
		
		TextWindow textWindow= new TextWindow();
		textWindow.setFrom(101);
		textWindow.setTo(120);

		SAMRecord rec= new SAMRecord(null);
		rec.setAlignmentStart(101); rec.setReadString("AAAAA"); rec.setCigarString("5M"); 
		SAMRecord rec2= new SAMRecord(null);
		rec2.setAlignmentStart(102); rec2.setReadString("CCCCC"); rec2.setCigarString("5M"); 
		SAMRecord rec3= new SAMRecord(null);
		rec3.setAlignmentStart(106); rec3.setReadString("TTTTT"); rec3.setCigarString("5M"); 
		SAMRecord rec4= new SAMRecord(null);
		rec4.setAlignmentStart(109); rec4.setReadString("GGGGG"); rec4.setCigarString("5M"); 
		SAMRecord rec5= new SAMRecord(null);
		rec5.setAlignmentStart(115); rec5.setReadString("AAAAA"); rec5.setCigarString("5M"); 
		
		List<TextRead> textReads= new ArrayList<TextRead>();
		textReads.add(new TextRead(rec, textWindow));
		textReads.add(new TextRead(rec2, textWindow));
		textReads.add(new TextRead(rec3, textWindow));
		textReads.add(new TextRead(rec4, textWindow));
		textReads.add(new TextRead(rec5, textWindow));

		System.out.println("Unstacked:");
		for(TextRead tr : textReads){
			System.out.println(tr.getPrintableText());
		}
		System.out.println("\nStacked:");
		List<List<TextRead>> stackReads= Utils.stackReads(textReads);
		for(List<TextRead> line : stackReads){
			System.out.println(Utils.linePrinter(line, true));
		} 
		
//		System.out.println("\nPrint one read as bold:");
//		rec.setReadNegativeStrandFlag(true);
//		TextRead tr= new TextRead(rec, textWindow);
//		List<TextRead> trList= new ArrayList<TextRead>();
//		trList.add(tr);
//		System.out.println(Utils.linePrinter(trList, true));
		
	}
	
	@Test 
	public void canGetStartOfBamFile(){
		GenomicCoords gc= Utils.getStartCoordsOfBAM("test/test_data/ds051.short.bam");
		assertEquals("chr7", gc.getChrom());
		assertEquals(5566778, (int)gc.getFrom());
		assertEquals(gc.getFrom(), gc.getTo());

	}
	
	@Test
	public void canGetStartOfChrom(){
		GenomicCoords gc= Utils.getStartCoordsOfBAM("test/test_data/ds051.short.bam", "chr7");
		assertEquals("chr7", gc.getChrom());
		assertEquals(5566778, (int)gc.getFrom());
		assertEquals(gc.getFrom(), gc.getTo());
		
		// null for chrom with no reads. Fails if chrom does not exist in header.
		gc= Utils.getStartCoordsOfBAM("test/test_data/ds051.short.bam", "chrY");
		assertEquals(null, gc.getChrom());
		assertEquals(null, gc.getFrom());
		assertEquals(null, gc.getTo());
	}
	
	@Test
	public void canParseInputAndUpdateGenomicCoords(){
		GenomicCoords gc= new GenomicCoords("chr7:100-200", samSeqDict);

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
	
	@Test
	public void canCompressListOfInts(){
		
		// N. windows is a multiple of N. elements
		int[] intArr= {1, 2, 3, 10, 20, 30, 100, 200, 300};
		List<Integer> intList= new ArrayList<Integer>();
		for(int x : intArr){
			intList.add(x);
		}
		int nwinds= 3;
		LinkedHashMap<Integer, Integer> zlist= Utils.compressListOfInts(intList, nwinds);
		Set<Integer> at= zlist.keySet();
		assertEquals("[0, 3, 6]", at.toString());
		Collection<Integer> depth= zlist.values();
		assertEquals("[2, 20, 200]", depth.toString());

		// No compression as more windows then elements
		nwinds= 300; 
		zlist= Utils.compressListOfInts(intList, nwinds);
		assertEquals(intList, new ArrayList<Integer>(zlist.values()));
		
		// An odd division of #elements by #windows
		intList.clear();
		for (int i = 0; i < 100; i++) {
			intList.add(i);
		}
		nwinds= 17;
		zlist= Utils.compressListOfInts(intList, nwinds);
		assertEquals(nwinds, zlist.size());
		assertEquals("[3, 9, 15, 21, 27, 33, 39, 45, 51, 57, 63, 69, 75, 81, 87, 93, 98]", zlist.values().toString());

		// Empty input list return empty map. 
		nwinds= 3; 
		intList.clear();
		zlist= Utils.compressListOfInts(intList, nwinds);
		assertEquals(0, zlist.size());
	}
	
	@Test
	public void canPrintRuler(){
		int from= 1;
		int to= 1000;
		int by= 10;
		int windowSize= 160;
		String ruler= Utils.ruler(from, to, by, windowSize);
		System.out.println("Ruler:");
		System.out.println(ruler);
	}
}