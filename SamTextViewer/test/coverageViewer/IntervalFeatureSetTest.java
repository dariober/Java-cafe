package coverageViewer;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.lang.instrument.Instrumentation;
import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import com.google.common.collect.ContiguousSet;
import com.google.common.collect.DiscreteDomain;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Range;

public class IntervalFeatureSetTest {

	@Test
	public void test() throws IOException {
		// Note that refSeq.hg19.short.bed is not sorted by pos. 
		IntervalFeatureSet set= new IntervalFeatureSet(new File("test_data/refSeq.hg19.short.bed"));
		assertEquals(set.getIntervalMap().get("chr1").get(0).getName(), "NM_001080397_utr3_8_0_chr1_8404074_f");
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
	
	@Test
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
		List<Float> rulerMap= new ArrayList<Float>();
		for(int i= 30000000; i < 70000000; i += 2000000){
			rulerMap.add((float)(i + 0.3));
		} // [3.0E7, 3.2E7, 3.4E7, 3.6E7, 3.8E7, 4.0E7, 4.2E7, 4.4E7, 4.6E7, 4.8E7, 5.0E7, 5.2E7, 5.4E7, 5.6E7, 5.8E7, 6.0E7, 6.2E7, 6.4E7, 6.6E7, 6.8E7]
		IntervalFeatureSet set= new IntervalFeatureSet(new File("test_data/refSeq.hg19.short.bed"));
		set.mapIntervalsToScreen("chr1", rulerMap);
		
		List<IntervalFeature> mapSet = set.getIntervalMap().get("chr1");
		
		assertEquals(-1, mapSet.get(0).getScreenFrom());
		assertEquals(19, mapSet.get(mapSet.size() - 1).getScreenFrom());
		assertEquals(19, mapSet.get(mapSet.size() - 1).getScreenTo());
	}	
	
	@Test
	public void canPrintMappingOfFeaturesToScreen() throws IOException{
		List<Float> rulerMap= new ArrayList<Float>();
		for(int i= 14000; i < 14400; i += 10){
			rulerMap.add((float)i);
		}
		System.out.println(rulerMap);
		IntervalFeatureSet set= new IntervalFeatureSet(new File("test_data/refSeq.hg19.bed.gz"));
		System.out.println(set.getIntervalMap().get("chr1").get(0));
		System.out.println(set.printableIntervals("chr1", rulerMap));
	}
}
