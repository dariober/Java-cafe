package samTextViewer;

import org.apache.commons.lang3.StringUtils;

import readWriteBAMUtils.ReadWriteBAMUtils;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.ValidationStringency;

/**Parse and store genomic coordinates 
 * TODO: Commas from ints are stripped, regardless of whether they are in the right place or not. 
 * @author berald01
 */
public class GenomicCoords {
	
	private String chrom;
	private Integer from;
	private Integer to;
	private SAMSequenceDictionary samSeqDict;
	
	/* Constructors */
	public GenomicCoords(String region, SAMSequenceDictionary samSeqDict){
		
		parseStringToGenomicCoords(region);
		this.samSeqDict= samSeqDict;
		makeGenomicCoordsValid();
		correctCoordsAgainstSeqDict(samSeqDict);

	}
	
	public GenomicCoords(String chrom, Integer from, Integer to, SAMSequenceDictionary samSeqDict){
		this.samSeqDict= samSeqDict;
		this.chrom= chrom;
		this.from= from;
		this.to= to;
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
		return this.chrom + ":" + this.from + "-" + this.to;
	}
	
	/* Getters and setters */
	
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