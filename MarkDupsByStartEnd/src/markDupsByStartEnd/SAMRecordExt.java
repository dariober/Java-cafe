package markDupsByStartEnd;

import htsjdk.samtools.SAMRecord;

/**
 * @author berald01
 * Extended SAMRecord class to be used to appropriately sort records for
 * duplicate marking. Hold basequality score, RG field, unclipped start and end etc...  
 */
public class SAMRecordExt implements Comparable<SAMRecordExt>{
	
	private short baseQualityScore= -1;
	private SAMRecord samrecord;
	private String rgtag;
	//private Integer referenceIndex;
	//private int unclippedStart;
	//private int unclippedEnd;
	//private boolean negStrand;
	//private int mapq;
	
	/* C O N S T R U C T O R */
	public SAMRecordExt() {

	}	

	public SAMRecordExt(SAMRecord rec, boolean ignoreReadGroup) {
		
		//this.referenceIndex= rec.getReferenceIndex();
		//this.unclippedStart= rec.getUnclippedStart();
		//this.unclippedEnd= rec.getUnclippedEnd();
		//this.negStrand= rec.getReadNegativeStrandFlag();
		this.rgtag= getRGtag(rec, ignoreReadGroup);
		this.baseQualityScore= getSumOfBaseQualities(rec);
		//this.mapq= rec.getMappingQuality();
		this.samrecord= rec;
		
	}
	
	public void setBaseQualityScore(short baseQualityScore) {
		this.baseQualityScore = baseQualityScore;
	}
	
	/** Calculates a score for the read which is the sum of scores over Q15. 
	 * Code taken from 
	 * https://github.com/samtools/htsjdk/blob/master/src/java/htsjdk/samtools/DuplicateScoringStrategy.java
	 * */
	public static short getSumOfBaseQualities(SAMRecord samrecord) {
		short score = 0;
		for (byte b : samrecord.getBaseQualities()) {
			if (b >= 15) score += b;
			}
		return score;
	}
	
	/**
	 * Return a string of read group ID or an not-available character if RG is
	 * not avalilable or not required.
	 * @param samrecord
	 * @param ignoreReadGroup Should the read group be ignored? If true an NA string is returned.
	 * @return
	 */
	private static String getRGtag(SAMRecord samrecord, boolean ignoreReadGroup){
		// Use '\b' as NA since it is not allowed in RG names.
		String rgtag= (samrecord.getAttribute("RG") == null || ignoreReadGroup) ? "\b" : samrecord.getReadGroup().getId();
		return rgtag;
	}

	public String toString(){
		// TODO: Return object as string suitable to be read line by line later.
		return this.toString();
	}
	
	/*  Comparator  */	
	public int compareTo(SAMRecordExt other) {
	    int i = this.samrecord.getReferenceIndex()- other.samrecord.getReferenceIndex();
	    if (i != 0) return i;
	    
	    i = this.samrecord.getUnclippedStart() - other.samrecord.getUnclippedStart();
	    if (i != 0) return i;
	    
	    i = this.samrecord.getUnclippedStart() - other.samrecord.getUnclippedStart();
	    if (i != 0) return i;
	    
	    i= Boolean.valueOf(this.samrecord.getReadNegativeStrandFlag()).compareTo(other.samrecord.getReadNegativeStrandFlag());
	    if (i != 0) return i;
	    
	    i= this.rgtag.compareTo(other.rgtag);
	    if (i != 0) return i;
	    
	    i= other.baseQualityScore - this.baseQualityScore; // DESC: other - this
	    if (i != 0) return i;
	    
	    i= other.samrecord.getMappingQuality() - this.samrecord.getMappingQuality();
	    if (i != 0) return i;
		
	    return i;
	    
	}
	
	/* S E T T E R S  A N D  G E T T E R S */
	
	public SAMRecord getSamrecord() {
		return samrecord;
	}

	public void setSamrecord(SAMRecord samrecord) {
		this.samrecord = samrecord;
	}

	public String getRgtag() {
		return rgtag;
	}

	public void setRgtag(String rgtag) {
		this.rgtag = rgtag;
	}

	public short getBaseQualityScore() {
		return baseQualityScore;
	}
	
}
