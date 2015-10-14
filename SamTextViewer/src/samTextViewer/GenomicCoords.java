package samTextViewer;

import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;

import org.apache.commons.lang3.StringUtils;

import readWriteBAMUtils.ReadWriteBAMUtils;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;

/**
 * Class to set up the horizontal axis on screen. 
 * 
 * * Parse and check input genomic coordinates in input.
 * * Map genomic coordinates to screen coordinates.
 * * Display info about window: Ruler_TO_BE_DEPRECTED, coordinates, bp/column
 * 
 * TODO: Do not check chrom and upper limit if a SamFileHeader is not given, e.g.
 * because ot was not loaded.
 *  
 * @author berald01
 */
public class GenomicCoords {
	
	private String chrom;
	private Integer from;
	private Integer to;
	private SAMSequenceDictionary samSeqDict;
	// private List<Double> mapping= new ArrayList<Double>();
	
	/* Constructors */
	public GenomicCoords(String region, SAMSequenceDictionary samSeqDict){
		
		parseStringToGenomicCoords(region);
		this.samSeqDict= samSeqDict;
		makeGenomicCoordsValid();
		correctCoordsAgainstSeqDict(samSeqDict);

	}
	
	public GenomicCoords(String chrom, Integer from, Integer to, SAMSequenceDictionary samSeqDict){
		
		this.chrom= chrom;
		this.from= from;
		this.to= to;		
		this.samSeqDict= samSeqDict;
		makeGenomicCoordsValid();
		correctCoordsAgainstSeqDict(samSeqDict);

	}
	
	/* Methods */
	
	/**
	 * Parse string to return a GenomicCoords obj. This object
	 * can't be used as such as there is no check for valid input
	 * @return
	 */
	private void parseStringToGenomicCoords(String x){
		
		x= x.trim();
		Integer from= null;
		Integer to= null;
		int nsep= StringUtils.countMatches(x, ":");
		if(nsep == 0){ // Only chrom present
			this.chrom= x.trim();
			this.from= from;
			this.to= to;
			return;
		} else {

			this.chrom= StringUtils.substringBeforeLast(x, ":").trim();
		
			String fromTo= StringUtils.substringAfterLast(x, ":").replaceAll(",", "").trim();
			if(!fromToRegionIsValid(fromTo)){
				this.from= null;
				this.to= null;
				return;		
			}
			nsep= StringUtils.countMatches(fromTo, "-");
			if(nsep == 0){
				this.from= Integer.parseInt(StringUtils.substringBefore(fromTo, "-").trim());
				this.to= null;
				return;
			} else if(nsep == 1){
				this.from= Integer.parseInt(StringUtils.substringBefore(fromTo, "-").trim());
				this.to= Integer.parseInt(StringUtils.substringAfter(fromTo, "-").trim());
				return;
			} else {
				System.err.println("Unexpected format for region " + x + " returning null");				
				System.exit(1);
			}
		} 
	} 

	/**
	 * Make sure that the string after ":" in chr:start-end
	 * is valid. 
	 * @param x
	 * @return
	 */
	private boolean fromToRegionIsValid(String x){
		x= x.trim();
		if(x.length() == 0){
			return false;
		} else if(x.startsWith("-")){
			return false;
		} else if(StringUtils.countMatches(x, "-") > 1){
			return false;
		} else {
			x= x.replaceAll("\\s*\\-\\s*", "");
			for(int i=0; i < x.length(); i++){
				if( ! Character.isDigit(x.charAt(i))){
					return false;
				}
			}
			if(x.replaceAll("0", "") == ""){ // Region is not only zeros!
				return false;
			}
		}
		return true;
	}
	
	private void makeGenomicCoordsValid(){
		
		// Checks from chrom:
		if(this.chrom == null){ // Chrom must be not null otherwise nullify everything
			this.from= null;
			this.to= null;
			return;
		} 
		// Passed this point chrom can't be null.

		// Checks for interval:
		if(this.from == null){
			this.to= null;
			return;
		}

		if(this.from <= 0){
			this.from= null;
			this.to= null;
			return;			
		}
		
		if(this.to != null && this.from > this.to){
			this.from= null;
			this.to= null;
			return;						
		}
		
	}
	
	public void correctCoordsAgainstSeqDict(SAMSequenceDictionary samSeqDict){
		
		if(this.chrom == null){ // Nothing to do
			return;
		}
		if(samSeqDict.getSequence(this.chrom) == null){ // Not found: Nullify everything
			this.chrom= null;
			this.from= null;
			this.to= null;
			return;
		} 
		// Reset max coords
		if( this.from != null && this.from > samSeqDict.getSequence(this.chrom).getSequenceLength() ) {
			this.from= samSeqDict.getSequence(this.chrom).getSequenceLength();
		}
		if( this.to != null && this.to > samSeqDict.getSequence(this.chrom).getSequenceLength() ) {
			this.to= samSeqDict.getSequence(this.chrom).getSequenceLength();
		}
	}

	/** 
	 * Go to region specified by the given coordinates. Adjust coordinates as neccessary
	 * @param windowSize Reset coordinates to fit this window size if only chrom or chrom:start
	 * is not null.
	 * @return A new GenomicRegion object with coordinate moved to region. Adjusted 
	 * to be of size windowSize.
	 */
	public static GenomicCoords goToRegion(String region, String bam, int windowSize){
		
		SAMSequenceDictionary samSeqDict= ReadWriteBAMUtils
				.reader(bam, ValidationStringency.LENIENT)
				.getFileHeader()
				.getSequenceDictionary();
		
		GenomicCoords gc= new GenomicCoords(region, samSeqDict);
		
		Integer chromSize= null;
		if(gc.chrom != null){
			chromSize= samSeqDict.getSequence(gc.chrom).getSequenceLength();
		} else {
			gc.correctCoordsAgainstSeqDict(samSeqDict);
			return gc;
		}
		
		if(gc.from != null && gc.to != null){
			
			if((gc.from + windowSize - 1) > chromSize){
				gc.from= chromSize - windowSize + 1;
				gc.to= gc.from + windowSize - 1;
			}			
		} else if(gc.from != null && gc.to == null){
			if((gc.from + windowSize - 1) > chromSize){
				gc.from= chromSize - windowSize + 1;
			}
			gc.to= gc.from + windowSize - 1;
		} else if(gc.from == null){
			GenomicCoords gcStart= Utils.getStartCoordsOfBAM(bam, gc.chrom);
			if(gcStart.getFrom() == null){ // This is if there are no mapped reads
				gc.from= 1;
				gc.to= windowSize;
			} else {
				gc.from= gcStart.getFrom();
				gc.to= gcStart.getFrom() + windowSize - 1;
			}
		}
		gc.correctCoordsAgainstSeqDict(samSeqDict);
		return gc;
	}
	
	public String toString(){
		int range= this.to - this.from + 1;
		return this.chrom + ":" + this.from + "-" + this.to + "; " + NumberFormat.getNumberInstance(Locale.UK).format(range) + " bp";
	}
	
	private int getMidpoint(){
		int range= this.to - this.from + 1;
		if(range % 2 == 1){
			range--;
		}
		// * Get midpoint of genomic interval
		int midpoint= range / 2 + this.from;
		return midpoint;
	}
	
	/**
	 * Rescale coords to extend them as in zooming-in/-out
	 */
	public void zoomOut(){
		int zoom= 2;
		// * Get size of window (to - from + 1)
		int range= this.to - this.from + 1;
		if(range % 2 == 1){
			range--;
		}
		int midpoint= this.getMidpoint();
		
		// * Extend midpoint left by window size x2 and check coords
		this.from= midpoint - (range * zoom);
		this.from= (this.from <= 0) ? 1 : this.from; 
		
		// Extend midpoint right
		this.to= midpoint + (range * zoom);
		if(this.chrom != null && this.samSeqDict.getSequence(this.chrom).getSequenceLength() > 0){
			this.to= (this.to > this.samSeqDict.getSequence(this.chrom).getSequenceLength()) ? 
					this.samSeqDict.getSequence(this.chrom).getSequenceLength() : this.to;
		}
	}

	/**
	 * Zoom into range. 
	 */
	public void zoomIn(){
		float zoom= (float) (1/4.0);
		// * Get size of window (to - from + 1)
		int range= this.to - this.from + 1;
		if(range % 2 == 1){
			range--;
		}
		// * Get midpoint of range
		int midpoint= this.getMidpoint();
		int extendBy= Math.round(range * zoom);
		this.from= midpoint - extendBy;
		this.to= midpoint + extendBy;
		if(this.from > this.to){ // Not sure this can happen.
			this.to= this.from;
		}
	}
	
	/**
	 * Same as R seq(from to, length.out). 
	 * @param from
	 * @param to
	 * @param lengthOut Length of sequence, effectively the desired screen width.
	 * @return
	 */
	private List<Double> seqFromToLenOut(int lengthOut){
		
		List<Double> mapping= new ArrayList<Double>();
		
		if(from < 1 || from > to){
			System.err.println("Invalid genome coordinates: from " + from + " to " + to);
			System.exit(1);
		}
		int span= to - from + 1;
		// If the genomic span is less then screen size, reduce screen size to.
		// If genomic span == screenSize then you have a mapping one to one.
		if(span <= lengthOut){ 
			for(int i= from; i <= to; i++){
				mapping.add((double)i);
			}
			return mapping;
		}
		
		double step= ((double)span - 1)/(lengthOut - 1);
		mapping.add((double)from);
		for(int i= 1; i < lengthOut; i++){
			mapping.add((double)mapping.get(i-1)+step);
		}
		
		// First check last point is close enough to expection. If so, replace last point with
		// exact desired.
		double diffTo= Math.abs(mapping.get(mapping.size() - 1) - to);
		if(diffTo > ((float)to * 0.001)){
			System.err.println("Error generating sequence:");
			System.err.println("Last point: " + mapping.get(mapping.size() - 1));
			System.err.println("To diff: " + diffTo);
			System.err.println("Step: " + step);
		} else {
			mapping.set(mapping.size()-1, (double)to);
		}
		
		double diffFrom= Math.abs(mapping.get(0) - from);		
		if(diffFrom > 0.01 || mapping.size() != lengthOut){
			System.err.println("Error generating sequence:");
			System.err.println("Expected size: " + lengthOut + "; Effective: " + mapping.size());
			System.err.println("From diff: " + diffFrom);
			System.exit(1);
		}
		return mapping;
	}
	
	public double getBpPerScreenColumn(int windowSize){
		List<Double> mapping = seqFromToLenOut(windowSize);
		double bpPerScreenColumn= (to - from + 1) / (double)mapping.size();
		return bpPerScreenColumn;
	}
	
	/**
	 * Mapping of genome positions to screen. Screen position is
	 * 0-based. I.e. first column has position 0. For genomic
	 * spans only slightly larger then screen size the mapping might be
	 * inaccurate due to rounding during preparation of windows. 
	 * */
	public int getScreenPositionAtGenomePosition(int windowSize, int genomePos){
		List<Double> mapping = seqFromToLenOut(windowSize);
		if(genomePos < from || genomePos > to){
			return -1;
		}
		int closest= Utils.getIndexOfclosestValue((double)genomePos, mapping);
		return closest;
	}
	
	/** For debugging only */
	public String toStringVerbose(int windowSize){
		List<Double> mapping = seqFromToLenOut(windowSize);
		String str= "Genome coords: " + from + "-" + to 
				+ "; screen width: " + mapping.size()
				+ "; scale: " + this.getBpPerScreenColumn(windowSize) + " bp/column" 
				+ "; Mapping: " + mapping;
		str += "\n";
		str += this.toString();
		return str;
	}
	
	public String printableRuler(int windowSize, int markDist){
		List<Double> mapping = seqFromToLenOut(windowSize);
    	String numberLine= "";
    	int prevLen= 0;
    	int i= 0;
		while(i < mapping.size()){
			String posMark= String.valueOf(Math.round(mapping.get(i)));
			if(i == 0){
				numberLine= posMark;
				i += posMark.length();
			} else if((numberLine.length() - prevLen) >= markDist){
				prevLen= numberLine.length();
				numberLine= numberLine + posMark;
				i += posMark.length();
			} else {
				numberLine= numberLine + " ";
				i++;
			}
		}
    	return numberLine;	
    }
	
	/* Getters and setters */
	
	public List<Double> getMapping(int windowSize) {
		return seqFromToLenOut(windowSize);
	}	
	
	public String getChrom() {
		return chrom;
	}

	public void setChrom(String chrom) {
		this.chrom = chrom;
	}

	public Integer getFrom() {
		return from;
	}

	public void setFrom(int from) {
		if(from <= 0){
			from= 1;
		}
		this.from = from;
	}

	public Integer getTo() {
		return to;
	}

	public void setTo(int to) {
		if(to <= 0){
			to= 1;
		}
		this.to = to;
	}

	public SAMSequenceDictionary getSamSeqDict() {
		return samSeqDict;
	}

	public void setSamSeqDict(SAMSequenceDictionary samSeqDict) {
		this.samSeqDict = samSeqDict;
	}
	
}