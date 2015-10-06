package coverageViewer;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Random;

import org.apache.commons.lang3.StringUtils;
import org.junit.Test;

import samTextViewer.Utils;

public class CoverageViewerTest {

	@Test
	public void stubCovergeViewer(){
		
		String sam= "test_data/ds051.short.bam";
		String chrom= "chr7";
		int from= 5566781;
		int to= from+10;
		
		CoverageViewer cw= new CoverageViewer(sam, chrom, from, to, 1000);
		assertEquals(to - from + 1, cw.getDepth().size());
		assertEquals("[4, 5, 6, 8, 9, 9, 9, 9, 9, 9, 10]", cw.getDepth().toString());
		assertEquals(10, cw.getMaxDepth());
		assertEquals(5566781, (int)cw.getDepthAt().get(0));
		
		String expProfile= ""
				+ "          *\n"
				+ "    *******\n"
				+ "   ********\n"
				+ "   ********\n"
				+ "  *********\n"
				+ " **********\n"
				+ "***********\n"
				+ "***********\n"
				+ "***********\n"
				+ "***********";
		//List<String> depthList = 
		
		String obsProfile= StringUtils.join(cw.getProfileStrings(10), "\n");
		assertEquals(expProfile, obsProfile);

		expProfile= ""
				+ "    *******\n"
				+ "   ********\n"
				+ " **********\n"
				+ "***********\n"
				+ "***********";		
		obsProfile= StringUtils.join(cw.getProfileStrings(5), "\n");
		assertEquals(expProfile, obsProfile);
		
	}
	
	@Test
	public void canCompressCoverageView(){
		
		String sam= "test_data/ds051.short.bam";
		String chrom= "chr7";
		int from= 5566781;
		int to= from+200;
		int windowSize= 10;
		
		CoverageViewer cw= new CoverageViewer(sam, chrom, from, to);
		cw.compressCovergeViewer(windowSize);
		cw.getProfileStrings(10);
		String expRuler= "5566781---5566981";
		assertEquals(expRuler, cw.ruler(10));		
	}
	
	@Test
	public void canConstructLargeIntervalBySampling(){
		
		String sam= "test_data/mjb050_oxBS.bam";
		String chrom= "chrY";
		int from= 10000;
		int to= 50000;
		int windowSize= 50;
		CoverageViewer cw= new CoverageViewer(sam, chrom, from, to, windowSize);
		// The number of loci collected is ~ LOC_PER_WINDOW * windowSize
		assertEquals(CoverageViewer.LOC_PER_WINDOW * windowSize, cw.getDepth().size(), 200);
		
		// Start sampling at approx start of coords
		assertEquals(from, cw.getDepthAt().get(0), 50);
		// Stop sampling at approx end of coords.
		assertEquals(to, cw.getDepthAt().get(cw.getDepthAt().size()-1), 50);
	}
	
	@Test
	public void canCompressLargeCoverageView(){
		
		String sam= "test_data/ds051.actb.bam";
		String chrom= "chr7";
		int from= 5566781;
		int to= from+10000;
		int windowSize= 50;
		
		CoverageViewer cw= new CoverageViewer(sam, chrom, from, to, windowSize);
		cw.compressCovergeViewer(windowSize);
		assertEquals(from, cw.getDepthAt().get(0), 200);
		assertEquals(to, cw.getDepthAt().get(cw.getDepthAt().size()-1), 200);
//		System.out.println(cw);
	}

}

/*
 		long t0= System.currentTimeMillis();
		long t1= System.currentTimeMillis();
		System.out.println(t1-t0 + "ms");
*/
