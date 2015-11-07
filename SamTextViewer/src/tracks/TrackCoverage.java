package tracks;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.apache.commons.lang3.StringUtils;

import com.google.common.base.Joiner;

import htsjdk.samtools.BAMIndex;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import samTextViewer.GenomicCoords;
import samTextViewer.SamLocusIterator;
import samTextViewer.Utils;

public class TrackCoverage {

	/* A t t r i b u t e s */
	
	/** Max number of loci per window (column of text char). If you query 
	 * large intervals, only some loci will be stored. This const define the density of
	 * loci per window. */
	// final static int LOC_PER_WINDOW= 200;
	
	
	private List<ScreenLocusInfo> screenLocusInfoList= new ArrayList<ScreenLocusInfo>(); 

	/** Each dot in the screen output track corresponds to this many units of 
	 * score in the input. Typically this "reads per dot". */
	private double scorePerDot;   
	private double maxDepth;
	private String inputBam;

	/** Store collected loci in the region of interest */
	//@Deprecated
	//final private ArrayList<LocusInfo> locusInfoList= new ArrayList<LocusInfo>();
	
	/* C o n s t r u c t o r */
	
	/**
	 * Construct coverage track from bam alignment in the provided interval. 
	 * Loci will be sampled according to the size of the interval and the size of the printable screen. 
	 * @param bam Input bam file
	 * @param gc Interval to sample positions from
	 * @param windowSize The size of the screen in number of characters.
	 * @param filters Record filters to apply to input sam records.
	 * @param bs Should loci be parsed also as BS-Seq data? 
	 * @throws IOException 
	 */
	public TrackCoverage(String bam, GenomicCoords gc,
			List<SamRecordFilter> filters, boolean bs) throws IOException{
		
		this.inputBam= bam;
		SamReaderFactory srf=SamReaderFactory.make();
		srf.validationStringency(ValidationStringency.SILENT);
		SamReader samReader= srf.open(new File(bam));
		
		// SAMFileHeader fh= samReader.getFileHeader();
		//IntervalList il= locusSampler(fh, gc);
		IntervalList il= new IntervalList(samReader.getFileHeader());
		il.add(new Interval(gc.getChrom(), gc.getFrom(), gc.getTo()));
		SamLocusIterator samLocIter= new SamLocusIterator(samReader, il, true);
		samLocIter.setSamFilters(filters);
		Iterator<samTextViewer.SamLocusIterator.LocusInfo> iter= samLocIter.iterator();
	
		for(int i= 0; i < gc.getMapping().size(); i++){
			screenLocusInfoList.add(new ScreenLocusInfo());	
		}
		
		while(iter.hasNext()){			
			samTextViewer.SamLocusIterator.LocusInfo locusInfo= iter.next();
			int screenPos= Utils.getIndexOfclosestValue(locusInfo.getPosition(), gc.getMapping());
			byte refBase= '\0';
			if(gc.getRefSeq() != null){
				refBase= gc.getRefSeq()[screenPos];
			}
			screenLocusInfoList.get(screenPos).increment(locusInfo, refBase, bs);
		}
		samLocIter.close();
		samReader.close();
	}
	
	/* M e t h o d s */

	/**
	 * Printable coverage track. The height of the track in lines is `yMaxLines`.
	 * @param screenToGenomeMap List of genomic positions corresponding to each column on screen.
	 * @param yMaxLines
	 * @param rpm Should read counts be normalized by library size as Reads Per Million
	 * @return HashMapwith with keys/values the printable characteristics of the track. 
	 */
	public String printToScreen(List<Double> screenToGenomeMap, int yMaxLines, boolean rpm){
		
		List<Double> yValues= new ArrayList<Double>();
		for(ScreenLocusInfo x : screenLocusInfoList){
			yValues.add(x.getMeanDepth());
		}
		
		TextProfile textProfile= new TextProfile(yValues, yMaxLines);

		this.maxDepth= textProfile.getMaxDepth();
		this.scorePerDot= textProfile.getScorePerDot();
		if(rpm){
			long libSize= getAlignedReadCount(new File(this.inputBam));
			this.maxDepth= this.maxDepth / libSize * 1000000;
			this.scorePerDot= this.scorePerDot / libSize * 1000000;
		}
		
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
	
	/**
	 * Sample loci in the GenomciCoords interval so that they are evenly spaced given
	 * the screen windowSize.
	 * @param fh
	 * @param gc
	 * @param windowSize
	 * @return
	 */
	//private IntervalList locusSampler(SAMFileHeader fh, GenomicCoords gc){
	//	int range= gc.getTo() - gc.getFrom() + 1;
	//	// If this ration is >1 then every single locus in the interval is sampled
	//	double density= ((double)gc.getUserWindowSize() * LOC_PER_WINDOW) / range;  
	//	
	//	IntervalList il= new IntervalList(fh);
	//	for(int i= gc.getFrom(); i <= gc.getTo(); i++){
	//		Random rand= new Random(); // Really you shouldn't use rnd to filter loci.
	//		float p= rand.nextFloat();
	//		if(p < density){
	//			Interval interval= new Interval(gc.getChrom(), i, i);
	//			il.add(interval);
	//		}
	//	}
	//	return il;
	//}
        
    private long getAlignedReadCount(File bam){

    	SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);
		BAMIndex sr= new SAMFileReader(bam).getIndex();
		long alnCount= 0; 
		int i= 0;
		while(true){
			try{
				alnCount += sr.getMetaData(i).getAlignedRecordCount();
			} catch(NullPointerException e){
				break;
			}
			i++;
		}
		sr.close();
		return alnCount;
    }
    
    /* S e t t e r s   and   G e t t e r s */
    
    /* This method makes sense calling only after having set the profile. Typically after  
     * */
    public double getScorePerDot(){
    	return this.scorePerDot;
    }
    public double getMaxDepth(){
    	return this.maxDepth;
    }    
    public List<ScreenLocusInfo> getScreenLocusInfoList(){
    	return screenLocusInfoList;
    }
    //public List<LocusInfo> getLocusInfoList(){
    //	return this.locusInfoList;
    //}
}

/**
 * Produce a representation of depth using text characters. 
 * Output is a list of lists where each inner list is a horizontal line of 
 * the coverage track. The vertical y-axis is scaled to be at most ymaxLines.
 * @param yValues Depth along positions
 * @param yMaxLines Max number of text lines to use. 
 * Values in depthArray will be rescaled accordingly. 0 or -ve to disable scaling.
 * @return
 */
/*protected List<List<String>> getProfileList(List<Double> yValues, int yMaxLines){
	
	double maxDepth= 0;
	for(double y : yValues){
		if(y > maxDepth){
			maxDepth= y;
		}
	}
	this.maxDepth= maxDepth;
	yMaxLines= yMaxLines * 2; // Since we use ':' for 2x in a single line.
	this.scorePerDot= (double)maxDepth/yMaxLines;
	
	// Rescale depth as required
	List<Double> yRescaled= new ArrayList<Double>();
	for(int i= 0; i < yValues.size(); i++){
		yRescaled.add(i, (double) yValues.get(i) / maxDepth * yMaxLines);  
	}
	
	List<List<String>> profile= new ArrayList<List<String>>();
	for(int i= 0; i < yRescaled.size(); i++){
		ArrayList<String> strDepth= new ArrayList<String>(); // This will be a vertical bar
		double locDepth= yRescaled.get(i);
		
		if((int)Math.rint(locDepth) == 0){ // For zero coverage
			strDepth.add("_"); 
		}
		int nDouble= ((int)locDepth) / 2; // how many :
		for(int j= 0; j < nDouble; j++){
			strDepth.add(":");
		} 
		int diff= (int)Math.rint(locDepth- (nDouble*2)); // Get remainder
		if(diff == 2){
			strDepth.add(":");
		} else if(diff == 1) {
			strDepth.add(".");
		} else if(diff == 0){
			//
		} else {
			System.err.println("Unexpected division");
			System.exit(1);
		}
		// Fill up list with blanks
		while(strDepth.size() < yMaxLines){
			strDepth.add(FILL);
	    }
		profile.add(strDepth);
	}
	return Utils.transpose(profile);
}*/

/**
 * Return map of positions and depth as extracted from the list of LocusInfo
 * */
/*private LinkedHashMap<Integer, Integer> getPositionsAndDepth(){   	
	LinkedHashMap<Integer, Integer> posAndDepth= new LinkedHashMap<Integer, Integer>();
	for(LocusInfo loc : this.locusInfoList){
		posAndDepth.put(loc.getPosition(), loc.getRecordAndPositions().size());
	}
	return posAndDepth;
}*/
