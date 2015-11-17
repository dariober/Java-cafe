package tracks;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
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
public class TrackWiggles extends Track {

	// private List<Double> screenScores= new ArrayList<Double>(); 
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

		this.setGc(gc);
		
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
			this.setScreenScores(screenScores);		
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
		for(int i= 0; i < getGc().getUserWindowSize(); i++){
			screenWiggleLocusInfoList.add(new ScreenWiggleLocusInfo());
		}
		
		BigWigIterator iter = reader.getBigWigIterator(getGc().getChrom(), getGc().getFrom(), getGc().getChrom(), getGc().getTo(), false);
		while(iter.hasNext()){
			WigItem bw = iter.next();
			for(int i= bw.getStartBase(); i <= bw.getEndBase(); i++){
				int idx= Utils.getIndexOfclosestValue(i, getGc().getMapping()); // Where should this position be mapped on screen?
				screenWiggleLocusInfoList.get(idx).increment(bw.getWigValue());
			} 
		}
		reader.close();
		ArrayList<Double> screenScores= new ArrayList<Double>();
		for(ScreenWiggleLocusInfo x : screenWiggleLocusInfoList){
			screenScores.add((double)x.getMeanScore());
		}
		this.setScreenScores(screenScores);		
	}
	
	private void bedGraphToScores(String fileName) throws IOException{
		
		List<ScreenWiggleLocusInfo> screenWiggleLocusInfoList= new ArrayList<ScreenWiggleLocusInfo>();
		for(int i= 0; i < getGc().getUserWindowSize(); i++){
			screenWiggleLocusInfoList.add(new ScreenWiggleLocusInfo());
		}
		
		// LineIterator it = FileUtils.lineIterator(new File(fileName), "UTF-8");
	
		try {
			TabixReader tabixReader= new TabixReader(fileName);
			Iterator qry= tabixReader.query(this.getGc().getChrom(), this.getGc().getFrom()-1, this.getGc().getTo());
			while(true){
				String q = qry.next();
				if(q == null){
					break;
				}
				String[] tokens= q.split("\t");
				int screenFrom= Utils.getIndexOfclosestValue(Integer.valueOf(tokens[1])+1, this.getGc().getMapping());
				int screenTo= Utils.getIndexOfclosestValue(Integer.valueOf(tokens[2]), this.getGc().getMapping());
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
		this.setScreenScores(screenScores);
		return;
	}
	
	@Override
	public String printToScreen(){ // int yMaxLines, Double ymin, Double ymax
		
		if(this.getyMaxLines() == 0){return "";}
		
		TextProfile textProfile= new TextProfile(this.getScreenScores(), this.getyMaxLines(), this.getYmin(), this.getYmax());
		
		this.scorePerDot= textProfile.getScorePerDot();
		this.maxDepth= textProfile.getMaxDepth();

		ArrayList<String> lineStrings= new ArrayList<String>();
		for(int i= (textProfile.getProfile().size() - 1); i >= 0; i--){
			List<String> xl= textProfile.getProfile().get(i);
			// Set<String> unique= new HashSet<String>(xl);
			lineStrings.add(StringUtils.join(xl, ""));
			//if(unique.size() == 1 && unique.contains(textProfile.getStrForFill())){ // Do not print blank lines
			//	continue;
			//} else {
			//	lineStrings.add(StringUtils.join(xl, ""));
			//}
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

	@Override
	public String getTitle(){
		return this.getFilename() + "; ylim: " + this.getYmin() + ", " + this.getYmax() + "; max: " + 
				Math.rint((this.maxDepth)*100)/100 + "; .= " + Math.rint((this.scorePerDot)) + ";\n";
	}
	
}
