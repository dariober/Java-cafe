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
			.reader("test_data/ds051.short.bam", ValidationStringency.STRICT)
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
		GenomicCoords gc= Utils.getStartCoordsOfBAM("test_data/ds051.short.bam");
		assertEquals("chr7", gc.getChrom());
		assertEquals(5566778, (int)gc.getFrom());
		assertEquals(gc.getFrom(), gc.getTo());

	}
	
	@Test
	public void canGetStartOfChrom(){
		GenomicCoords gc= Utils.getStartCoordsOfBAM("test_data/ds051.short.bam", "chr7");
		assertEquals("chr7", gc.getChrom());
		assertEquals(5566778, (int)gc.getFrom());
		assertEquals(gc.getFrom(), gc.getTo());
		
		// null for chrom with no reads. Fails if chrom does not exist in header.
		gc= Utils.getStartCoordsOfBAM("test_data/ds051.short.bam", "chrY");
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
		List<Float> longList= new ArrayList<Float>();
		for(float x : intArr){
			longList.add(x);
		}
		int nwinds= 3;
		LinkedHashMap<Integer, Float> compressedList= Utils.compressListOfInts(longList, nwinds);
		Set<Integer> at= compressedList.keySet();
		assertEquals("[0, 3, 6]", at.toString());
		Collection<Float> depth= compressedList.values();
		assertEquals("[2.0, 20.0, 200.0]", depth.toString());

		// No compression as more windows then elements
		nwinds= 300; 
		compressedList= Utils.compressListOfInts(longList, nwinds);
		assertEquals(longList, new ArrayList<Float>(compressedList.values()));
		
		// An odd division of #elements by #windows
		longList.clear();
		for (float i = 0; i < 100; i++) {
			longList.add(i);
		} // [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0, 60.0, 61.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0, 70.0, 71.0, 72.0, 73.0, 74.0, 75.0, 76.0, 77.0, 78.0, 79.0, 80.0, 81.0, 82.0, 83.0, 84.0, 85.0, 86.0, 87.0, 88.0, 89.0, 90.0, 91.0, 92.0, 93.0, 94.0, 95.0, 96.0, 97.0, 98.0, 99.0]

		System.out.println(longList);
		nwinds= 17;
		// System.out.println("Long list");
		compressedList= Utils.compressListOfInts(longList, nwinds);
		// System.out.println("END Long list");
		assertEquals(nwinds, compressedList.size());
		assertEquals("[2.5, 8.5, 14.5, 20.5, 26.5, 32.5, 38.5, 44.5, 50.5, 56.5, 62.5, 68.5, 74.5, 80.5, 86.5, 92.5, 98.0]", 
				compressedList.values().toString());
                    
		// BUG! compressedList.size() should be == nwinds but is rounded to longList.size()!
		// See also
		// echo -e "chr1\t0\t100" | windowMaker -b - -n 70 | awk '{print $0, $3-$2, NR}'
		nwinds= 70;
		compressedList= Utils.compressListOfInts(longList, nwinds);
		assertEquals(compressedList.size(), longList.size()); // 
		
		// Empty input list return empty map. 
		nwinds= 3; 
		longList.clear();
		compressedList= Utils.compressListOfInts(longList, nwinds);
		assertEquals(0, compressedList.size());
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