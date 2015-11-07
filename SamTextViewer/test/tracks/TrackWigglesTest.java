package tracks;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.List;

import org.broad.igv.tdf.TDFUtils;
import org.junit.Test;

import exceptions.InvalidGenomicCoordsException;
import htsjdk.tribble.readers.TabixReader;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;

public class TrackWigglesTest {

	@Test 
	public void canPrintBedGraph() throws InvalidGenomicCoordsException, IOException{
		
		String url= "test_data/test.bedGraph.gz";
		int yMaxLines= 5;
		int windowSize= 22;
		GenomicCoords gc= new GenomicCoords("chr1", 1, 22, null, windowSize, null);
		TrackWiggles tw= new TrackWiggles(url, gc);
		String prof= tw.printToScreen(yMaxLines);
		System.out.println(prof);
		assertEquals(5, tw.getMaxDepth(), 0.01);
		assertEquals(1, tw.getScorePerDot(), 0.01);
		
		tw= new TrackWiggles("test_data/positive.bedGraph.gz", gc);
		prof= tw.printToScreen(5);
		System.out.println(prof);
		assertEquals(5, tw.getMaxDepth(), 0.01);
		assertEquals(0.5, tw.getScorePerDot(), 0.01);
		
		tw= new TrackWiggles("test_data/negative.bedGraph.gz", gc);
		prof= tw.printToScreen(5);
		System.out.println(prof);
		assertEquals(0, tw.getMaxDepth(), 0.01);
		assertEquals(0.5, tw.getScorePerDot(), 0.01);
		
		gc= new GenomicCoords("chr1", 1, 52, null, 52, null);
		tw= new TrackWiggles("test_data/posNeg.bedGraph.gz", gc);
		System.out.println(tw.printToScreen(14));
	}
	

	
	@Test
	public void canPrintWiggleTrack() throws InvalidGenomicCoordsException, IOException {
		
		// * Check big* are 0 or 1 based
		
		// String url= "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCaltechRnaSeq/wgEncodeCaltechRnaSeqGm12892R2x75Th1014Il200SigRep3V4.bigWig";
		String url= "/Users/berald01/Downloads/wgEncodeCaltechRnaSeqGm12878R2x75Il400SigRep2V2.bigWig";
		
		int yMaxLines= 11;
		int windowSize= 101;
		GenomicCoords gc= new GenomicCoords("chrM", 1, 1000, null, windowSize, null);
		
		TrackWiggles tw= new TrackWiggles(url, gc);
		// System.out.println(tw.printToScreen(yMaxLines));
		
		//System.out.println(tw.getMaxDepth());
		//System.out.println(tw.getScorePerDot());
	
	}
	
	// @Test
	public void canPrintFromTdf() throws IOException, InvalidGenomicCoordsException{

		GenomicCoords gc= new GenomicCoords("chr8", 1, 100, null, 100, null);
		String tdfFile= "test_data/hg18_var_sample.wig.v2.1.30.tdf";
		List<ScreenWiggleLocusInfo> screenLocInfo = 
		TDFUtils.tdfRangeToScreen(tdfFile, gc.getChrom(), gc.getFrom(), gc.getTo(), gc.getMapping());
		// assertEquals(0.925, screenLocInfo.get(1).getMeanScore(), 0.1);

	
		gc= new GenomicCoords("chrM:1-16000", null, 100, null);
		tdfFile= "/Volumes/My_Passport_for_Mac/tmp/rhh_hacat_0508-1406_FAIRE.tdf";
		screenLocInfo = TDFUtils.tdfRangeToScreen(tdfFile, gc.getChrom(), gc.getFrom(), gc.getTo(), gc.getMapping());
		int i= 1;
		for(ScreenWiggleLocusInfo x : screenLocInfo){
			//System.out.println(i + " " + x);
			i++;
		}

		TrackWiggles tw= new TrackWiggles(tdfFile, gc);
		System.out.println(tw.printToScreen(40));
		System.out.println(tw.getMaxDepth());
		
	}

	
}
