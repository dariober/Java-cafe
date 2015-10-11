package coverageViewer;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.apache.commons.lang3.math.NumberUtils;
import org.junit.Test;

import samTextViewer.Utils;

public class IntervalFeatureTest {

	@Test
	public void beahviourOfIsNumber() {

		// Valid numbers
		assertTrue(NumberUtils.isNumber("1.1"));
		assertTrue(NumberUtils.isNumber("001"));
		assertTrue(NumberUtils.isNumber("0x0004")); // Also valid
		assertTrue(NumberUtils.isNumber("-1.1e9"));
		// Invalid numbers
		assertFalse(NumberUtils.isNumber("1.1 ")); // Note trailing space
		assertFalse(NumberUtils.isNumber(""));
		assertFalse(NumberUtils.isNumber("001.1"));
	}

	@Test
	public void canCreateIntervalFromString() {
		String bedLine= "chr1\t0\t1";
		IntervalFeature f= new IntervalFeature(bedLine, "bed");
		assertEquals("chr1", f.getChrom());
		assertEquals(1, f.getFrom()); // Note start augmented by 1.

		bedLine= "chr1\t0\t1\tgene\t0.1\t+";
		f= new IntervalFeature(bedLine, "bed");
		assertEquals("gene", f.getName());
		assertEquals('+', f.getStrand());

		bedLine= " chr1\t0\t1";
		f= new IntervalFeature(bedLine, "bed"); 
		assertEquals("chr1", f.getChrom()); // NB: spaces in chrom stripped.

		// System.out.println(f);
	}
	
	@Test
	public void canSortByChromPos() {
		List<IntervalFeature> flist= new ArrayList<IntervalFeature>();
		flist.add(new IntervalFeature("chrM\t10\t100\tg5", "bed"));
		flist.add(new IntervalFeature("chrM\t1\t100\tg4", "bed"));
		flist.add(new IntervalFeature("chrM\t1\t90\tg3", "bed"));
		flist.add(new IntervalFeature("chr1\t10\t90\tg2", "bed"));
		flist.add(new IntervalFeature("chr1\t1\t90\tg1", "bed"));
		Collections.sort(flist);
		assertEquals("g1", flist.get(0).getName());
		assertEquals("g5", flist.get(4).getName());
	}
	
	@Test
	public void canMapIntervalToRuler(){
		List<Float> rulerMap= new ArrayList<Float>();
		for(int i= 10; i < 20; i += 2){
			rulerMap.add((float)(i + 0.3));
		} // [10.3, 12.3, 14.3, 16.3, 18.3]
		IntervalFeature f= new IntervalFeature("chrM\t10\t15", "bed");
		f.mapToScreen(rulerMap);
		assertEquals(0, f.getScreenFrom());
		assertEquals(2, f.getScreenTo());

		// Feature is not part of ruler:
		f= new IntervalFeature("chrM\t100\t500", "bed");
		f.mapToScreen(rulerMap);
		assertEquals(-1, f.getScreenFrom());
		assertEquals(-1, f.getScreenTo());

		// NB: Feature is not part of ruler because start coord 18 is augmented by 1 to become 1-based.
		f= new IntervalFeature("chrM\t18\t30", "bed");
		f.mapToScreen(rulerMap);
		assertEquals(-1, f.getScreenFrom());
		assertEquals(-1, f.getScreenTo());
		
		// Partial overlap:
		f= new IntervalFeature("chrM\t1\t16", "bed");
		f.mapToScreen(rulerMap);
		assertEquals(0, f.getScreenFrom());
		assertEquals(3, f.getScreenTo());

		f= new IntervalFeature("chrM\t17\t30", "bed");
		f.mapToScreen(rulerMap);
		assertEquals(4, f.getScreenFrom());
		assertEquals(4, f.getScreenTo());
		//System.out.println(f);
		//System.out.println(rulerMap);
	}
}
