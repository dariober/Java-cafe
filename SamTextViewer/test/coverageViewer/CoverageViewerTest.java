package coverageViewer;

import static org.junit.Assert.*;

import java.io.File;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.lang3.StringUtils;
import org.junit.Test;

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
	
}
