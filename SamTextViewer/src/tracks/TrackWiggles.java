package tracks;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.lang3.StringUtils;
import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;
import org.broad.igv.tdf.TDFUtils;

import com.google.common.base.Joiner;

import htsjdk.tribble.readers.TabixReader;
import htsjdk.tribble.readers.TabixReader.Iterator;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;

/** Process wiggle file formats. Mostly using IGV classes. 
 * bigBed, bigWig, */
public class TrackWiggles {

	private List<Double> screenScores= new ArrayList<Double>(); 
	private GenomicCoords gc;
	// private double scorePerDot;   
	// private double maxDepth;
	private double maxDepth;
	private double scorePerDot;

	/* C o n s t r u c t o r s */

	/**
	 * Read bigWig from local file or remote URL.
	 * @param url Filename or URL to access 
	 * @param gc Query coordinates and size of printable window 
	 * @throws IOException */
	public TrackWiggles(String url, GenomicCoords gc) throws IOException{

		this.gc= gc;
		
		if(Utils.getFileTypeFromName(url).equals("bigWig")){
			BBFileReader reader=new BBFileReader(url); // or url for remote access.
			if(!reader.getBBFileHeader().isBigWig()){
				System.err.println("Invalid file type " + url);
				System.exit(1);			
			}
			bigWigToScores(reader);
		} else if(Utils.getFileTypeFromName(url).equals("tdf")){
			
			List<ScreenWiggleLocusInfo> screenWiggleLocusInfoList= TDFUtils.tdfRangeToScreen(url, gc.getChrom(), gc.getFrom(), gc.getTo(), gc.getMapping());
			ArrayList<Double> screenScores= new ArrayList<Double>();
			for(ScreenWiggleLocusInfo x : screenWiggleLocusInfoList){
				screenScores.add((double)x.getMeanScore());
			}
			this.screenScores= screenScores;		
		} else if(Utils.getFileTypeFromName(url).equals("bedGraph")){
			bedGraphToScores(url);
		} else {
			throw new RuntimeException("Extension (i.e. file type) not recognized for " + url);
		}
	};
	
	/*  M e t h o d s  */
	
	/** Populate object using bigWig data */
	private void bigWigToScores(BBFileReader reader){

		// List of length equal to screen size. Each inner map contains info about the screen locus 
		List<ScreenWiggleLocusInfo> screenWiggleLocusInfoList= new ArrayList<ScreenWiggleLocusInfo>();
		for(int i= 0; i < gc.getUserWindowSize(); i++){
			screenWiggleLocusInfoList.add(new ScreenWiggleLocusInfo());
		}
		
		BigWigIterator iter = reader.getBigWigIterator(gc.getChrom(), gc.getFrom(), gc.getChrom(), gc.getTo(), false);
		while(iter.hasNext()){
			WigItem bw = iter.next();
			for(int i= bw.getStartBase(); i <= bw.getEndBase(); i++){
				int idx= Utils.getIndexOfclosestValue(i, gc.getMapping()); // Where should this position be mapped on screen?
				screenWiggleLocusInfoList.get(idx).increment(bw.getWigValue());
			} 
		}
		reader.close();
		ArrayList<Double> screenScores= new ArrayList<Double>();
		for(ScreenWiggleLocusInfo x : screenWiggleLocusInfoList){
			screenScores.add((double)x.getMeanScore());
		}
		this.screenScores= screenScores;		
	}
	
	private void bedGraphToScores(String fileName){
		
		List<ScreenWiggleLocusInfo> screenWiggleLocusInfoList= new ArrayList<ScreenWiggleLocusInfo>();
		for(int i= 0; i < gc.getUserWindowSize(); i++){
			screenWiggleLocusInfoList.add(new ScreenWiggleLocusInfo());
		}
		try {
			TabixReader tabixReader= new TabixReader(fileName);
			Iterator qry= tabixReader.query(this.gc.getChrom(), this.gc.getFrom()-1, this.gc.getTo());
			while(true){
				String q = qry.next();
				if(q == null){
					break;
				}
				String[] tokens= q.split("\t");
				int screenFrom= Utils.getIndexOfclosestValue(Integer.valueOf(tokens[1])+1, this.gc.getMapping());
				int screenTo= Utils.getIndexOfclosestValue(Integer.valueOf(tokens[2]), this.gc.getMapping());
				float value= Float.valueOf(tokens[3]);
				for(int i= screenFrom; i <= screenTo; i++){
					screenWiggleLocusInfoList.get(i).increment(value);
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
			System.err.println("Could not open tabix file: " + fileName);
			System.err.println("Is the file sorted and indexed? After sorting by position (sort e.g. -k1,1 -k2,2n), compress with bgzip and index with e.g.:");
			System.err.println("\nbgzip " + fileName);
			System.err.println("tabix -p bed " + fileName + "\n");
		}
		ArrayList<Double> screenScores= new ArrayList<Double>();
		for(ScreenWiggleLocusInfo x : screenWiggleLocusInfoList){
			screenScores.add((double)x.getMeanScore());
		}
		this.screenScores= screenScores;
		return;
	}
	
	public String printToScreen(int yMaxLines){
		
		TextProfile textProfile= new TextProfile(this.screenScores, yMaxLines);
		
		this.scorePerDot= textProfile.getScorePerDot();
		this.maxDepth= textProfile.getMaxDepth();

		ArrayList<String> lineStrings= new ArrayList<String>();
		for(int i= (textProfile.getProfile().size() - 1); i >= 0; i--){
			List<String> xl= textProfile.getProfile().get(i);
			Set<String> unique= new HashSet<String>(xl);
			if(unique.size() == 1 && unique.contains(textProfile.getStrForFill())){ // Do not print blank lines
				continue;
			} else {
				lineStrings.add(StringUtils.join(xl, ""));
			}
		}
		return Joiner.on("\n").join(lineStrings);
	}
	
	/*   S e t t e r s   and   G e t t e r s */
	
	public double getMaxDepth() {
		return maxDepth;
	}

	public double getScorePerDot() {
		return scorePerDot;
	}


}
