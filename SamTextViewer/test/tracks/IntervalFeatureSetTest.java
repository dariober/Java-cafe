package tracks;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.lang.instrument.Instrumentation;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.apache.commons.lang3.text.StrTokenizer;
import org.junit.Test;

import com.google.common.collect.ContiguousSet;
import com.google.common.collect.DiscreteDomain;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Range;

import exceptions.InvalidGenomicCoordsException;
import samTextViewer.GenomicCoords;
import tracks.IntervalFeature;
import tracks.IntervalFeatureSet;

public class IntervalFeatureSetTest {

	@Test
	public void canReadFileWithHeader() throws IOException{
		IntervalFeatureSet set= new IntervalFeatureSet(new File("test_data/refSeq.bed"));
		assertEquals(2, set.getIntervalMap().size());
		assertEquals(2, set.getIntervalMap().get("chr1").size());
	}
	
	@Test
	public void test2() throws IOException{
		Set<String> chroms= new HashSet<String>();
		chroms.add("c1");
		chroms.add("c2");
		chroms.add("c3");
		chroms.add("c4");
		chroms.add("c5");
		chroms.add("c6");
		
		IntervalFeatureSet set= new IntervalFeatureSet(new File("test_data/refSeq.hg19.short.sort-2.bed"));
		System.out.println(set.getChromListStartingAt(chroms, "c3"));
		System.out.println(set.getChromListStartingAt(chroms, "c1"));
		System.out.println(set.getChromListStartingAt(chroms, "c6"));
	}

	
	@Test
	public void canFindNextFeatureOnChrom() throws IOException, InvalidGenomicCoordsException{

		//String tokens= "find  bla /my/  File";
		StrTokenizer str= new StrTokenizer("one\n   two '   three four'");
		str.setQuoteChar('\'');
		List<String> xs = str.getTokenList();		
		// new StringTokenizer("one two \"three four\"", ' ', '"').getTokenArray();
		
		System.out.println(xs);
		
		IntervalFeatureSet set= new IntervalFeatureSet(new File("test_data/refSeq.hg19.short.sort-2.bed"));
		
		IntervalFeature x = set.findNextStringOnChrom("nm_", "chr1", 20000000);
		assertTrue(x.getRaw().contains("NM_013943_utr3_5_0_chr1_25167429_f"));
		x = set.findNextStringOnChrom("nm_", "chr1", 80000000);
		assertTrue(x.getRaw().contains("NM_001080397_utr3_8_0_chr1_8404074_f"));
		
		x = set.findNextStringOnChrom("NotPresent", "chr1", 1);
		assertEquals(null, x);
		
		// Tabix
		set= new IntervalFeatureSet(new File("test_data/refSeq.hg19.short.sort.bed.gz"));
		x = set.findNextStringOnChrom("nm_", "chr1", 20000000);
		assertTrue(x.getRaw().contains("NM_013943_utr3_5_0_chr1_25167429_f"));
		x = set.findNextStringOnChrom("nm_", "chr1", 80000000);
		assertTrue(x.getRaw().contains("NM_001080397_utr3_8_0_chr1_8404074_f"));
	
		set= new IntervalFeatureSet(new File("test_data/refSeq.hg19.short.sort-2.bed"));
		x = set.findNextStringOnChrom("nm_", "chr1", 20000000);
		int i= 0;
		while(i < 20){
			// System.out.println(x);
			x = set.findNextStringOnChrom("nm_", x.getChrom(), x.getFrom());
			assertTrue(x.getRaw().contains("NM_"));
			i++;
		}

	}

	
	@Test
	public void canGetCoordsOfNextFeature() throws IOException, InvalidGenomicCoordsException{
		IntervalFeatureSet set= new IntervalFeatureSet(new File("test_data/refSeq.hg19.short.sort-2.bed"));
		
		GenomicCoords gc= new GenomicCoords("chr1:8000000-20000000", null, 100, null);
		GenomicCoords newGc= set.coordsOfNextFeature(gc);
		assertEquals(25167428+1, (int)newGc.getFrom());
		assertEquals(25167428+gc.getGenomicWindowSize(), (int)newGc.getTo());
	
		// Handling no next feature
		gc= new GenomicCoords("chr2:8000000-20000000", null, 100, null);
		newGc= set.coordsOfNextFeature(gc);
		assertEquals(gc, newGc);
		
		gc= new GenomicCoords("chr1:100000000-101000000", null, 100, null);
		newGc= set.coordsOfNextFeature(gc);
		assertEquals(gc, newGc);
		
	}

	@Test
	public void canGetCoordsOfNextFeatureTabix() throws IOException, InvalidGenomicCoordsException{
		IntervalFeatureSet set= new IntervalFeatureSet(new File("test_data/refSeq.hg19.short.sort.bed.gz"));
		
		GenomicCoords gc= new GenomicCoords("chr1:8000000-20000000", null, 100, null);
		GenomicCoords newGc= set.coordsOfNextFeature(gc);
		assertEquals(25167428+1, (int)newGc.getFrom());
		assertEquals(25167428+gc.getGenomicWindowSize(), (int)newGc.getTo());		

		gc= new GenomicCoords("chr2:8000000-20000000", null, 100, null);
		newGc= set.coordsOfNextFeature(gc);
		assertEquals(gc, newGc);
		
		gc= new GenomicCoords("chr1:100000000-101000000", null, 100, null);
		newGc= set.coordsOfNextFeature(gc);
		assertEquals(gc, newGc);
	}
	
	@Test
	public void test() throws IOException {
		// Note that refSeq.hg19.short.bed is not sorted by pos. 
		IntervalFeatureSet set= new IntervalFeatureSet(new File("test_data/refSeq.hg19.short.bed"));
		assertEquals("NM_001080397_utr3_8_0_chr1_8404074_f", set.getIntervalMap().get("chr1").get(0).getName());
		set= new IntervalFeatureSet(new File("test_data/refSeq.hg19.bed.gz"));		
	}

	@Test
	public void canFetchInterval() throws IOException{

		IntervalFeatureSet set= new IntervalFeatureSet(new File("test_data/refSeq.hg19.short.bed"));
		List<IntervalFeature> interval= set.getFeaturesInInterval("chr1", 20000000, 40000000);
		assertEquals(25167428+1, interval.get(0).getFrom()); // Note adding 1 because bed is 0-based
		assertEquals(33586132, interval.get(interval.size()-1).getTo());

		// Nothing to fetch: Range not in bed
		set= new IntervalFeatureSet(new File("test_data/refSeq.hg19.bed.gz"));
		interval= set.getFeaturesInInterval("chr1", 500000000, 600000000);
		assertEquals(0, interval.size());
		// Nothing to fetch: chrom not in bed:
		interval= set.getFeaturesInInterval("chrNonSense", 1, 10);
		assertEquals(0, interval.size());		
		
	}
	
	//@Test // BROKEN
	public void canMapIntervalSetToScreen() throws IOException{
		
/* Sorted refSeq.hg19.short.bed
chr1	8404073	8404227	NM_001080397_utr3_8_0_chr1_8404074_f	0	+
chr1	16785385	16786584	NM_001145278_utr3_7_0_chr1_16785386_f	0	+
chr1	16785385	16786584	NM_018090_utr3_7_0_chr1_16785386_f	0	+
chr1	16785491	16786584	NM_001145277_utr3_6_0_chr1_16785492_f	0	+
chr1	25167428	25170815	NM_013943_utr3_5_0_chr1_25167429_f	0	+
chr1	33585783	33586132	NM_001301825_utr3_8_0_chr1_33585784_f	0	+
chr1	33585783	33586132	NM_052998_utr3_11_0_chr1_33585784_f	0	+
chr1	33585783	33586132	NM_001293562_utr3_10_0_chr1_33585784_f	0	+
chr1	67208778	67216822	NM_032291_utr3_24_0_chr1_67208779_f	0	+
chr1	67208778	67216822	NM_001308203_utr3_21_0_chr1_67208779_f	0	+
*/
		List<Double> rulerMap= new ArrayList<Double>();
		for(int i= 30000000; i < 70000000; i += 2000000){
			rulerMap.add((double)(i + 0.3));
		} // [3.0E7, 3.2E7, 3.4E7, 3.6E7, 3.8E7, 4.0E7, 4.2E7, 4.4E7, 4.6E7, 4.8E7, 5.0E7, 5.2E7, 5.4E7, 5.6E7, 5.8E7, 6.0E7, 6.2E7, 6.4E7, 6.6E7, 6.8E7]
		IntervalFeatureSet set= new IntervalFeatureSet(new File("test_data/refSeq.hg19.short.bed"));
		//set.mapIntervalsToScreen("chr1", rulerMap);
		
		List<IntervalFeature> mapSet = set.getIntervalMap().get("chr1");
		
		assertEquals(-1, mapSet.get(0).getScreenFrom());
		assertEquals(19, mapSet.get(mapSet.size() - 1).getScreenFrom());
		assertEquals(19, mapSet.get(mapSet.size() - 1).getScreenTo());
	}	
	
	@Test
	public void canPrintMappingOfFeaturesToScreen() throws IOException{
		List<Double> rulerMap= new ArrayList<Double>();
		for(int i= 14000; i < 14400; i += 10){
			rulerMap.add((double)i);
		}
		System.out.println(rulerMap);
		IntervalFeatureSet set= new IntervalFeatureSet(new File("test_data/refSeq.hg19.bed.gz"));
		System.out.println(set.getIntervalMap().get("chr1").get(0));
	}
}
