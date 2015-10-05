package coverageViewer;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;

import org.apache.commons.lang3.StringUtils;
import org.junit.Test;

import samTextViewer.Utils;

public class CoverageViewerTest {

	@Test
	public void stubCovergeViewer(){
		
		String sam= "test/test_data/ds051.short.bam";
		String chrom= "chr7";
		int from= 5566781;
		int to= from+10;
		
		CoverageViewer cw= new CoverageViewer(sam, chrom, from, to);
		assertEquals(to - from + 1, cw.getDepth().size());
		assertEquals("[4, 5, 6, 8, 9, 9, 9, 9, 9, 9, 10]", cw.getDepth().toString());
		assertEquals(10, cw.getMaxDepth());
		assertEquals(5566781, (int)cw.getDepthAt().get(0));

//		cw.setDepthAt(cw.getDepthAt(), 1000000);
//		System.out.println(cw.getDepthAt());
		
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
		
		String sam= "test/test_data/ds051.short.bam";
		String chrom= "chr7";
		int from= 5566781;
		int to= from+200;
		int windowSize= 10;
		
		CoverageViewer cw= new CoverageViewer(sam, chrom, from, to);
		LinkedHashMap<Integer, Integer> zcw= Utils.compressListOfInts(cw.getDepth(), windowSize);
		ArrayList<Integer> zdepth= new ArrayList<Integer>(zcw.values());
		ArrayList<Integer> zat= new ArrayList<Integer>(zcw.keySet());
		cw= new CoverageViewer(zdepth, zat);
		cw.setDepthAt(cw.getDepthAt(), from);
		cw.getProfileStrings(10);
		System.out.println(cw);
		String expRuler= "5566781   5566801   5566821   5566841   5566861   5566881   5566901   5566921   5566941   5566961   5566981";
		System.out.println(cw.ruler(10));
//		cw.setDepthAt(cw.getDepthAt(), from);
		
	}
	
}
