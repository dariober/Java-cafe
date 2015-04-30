package markDupsByStartEnd;

import java.util.Arrays;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

/**
 * @author berald01
 * Extended SAMRecord class to be used to appropriately sort records for
 * duplicate marking. Hold basequality score, RG field, unclipped start and end etc...  
 */
public class SAMRecordExt implements Comparable<SAMRecordExt>{
	
	private short baseQualityScore= -1;
	private SAMRecord samRecord;
	private String rgtag;
	//private Integer referenceIndex;
	//private int unclippedStart;
	//private int unclippedEnd;
	//private boolean negStrand;
	//private int mapq;
	
	/* C O N S T R U C T O R */
	public SAMRecordExt() {

	}	

	/**
	 * Use a string, typically read from file, to produce a SAMRecord and 
	 * additional fields for SAMRecordExt.
	 * @param samRecordExtAsString String to be converted to record, tab separated.
	 * @param hdr
	 */
	public SAMRecordExt(String samRecordExtAsString, SAMFileHeader hdr) {
		String[] samRecordExtArray= samRecordExtAsString.split("\t");
		// Make sure these indexes are correct!!
		this.baseQualityScore= Short.parseShort(samRecordExtArray[0]); 
		this.rgtag= samRecordExtArray[1];
		String[] array= Arrays.copyOfRange(samRecordExtArray, 2, samRecordExtArray.length);
		this.samRecord= Utils.arrayToSAMRecord(array, hdr);
	}	
	
	public SAMRecordExt(SAMRecord rec, boolean ignoreReadGroup) {
		
		//this.referenceIndex= rec.getReferenceIndex();
		//this.unclippedStart= rec.getUnclippedStart();
		//this.unclippedEnd= rec.getUnclippedEnd();
		//this.negStrand= rec.getReadNegativeStrandFlag();
		this.rgtag= getRGtag(rec, ignoreReadGroup);
		this.baseQualityScore= getSumOfBaseQualities(rec);
		//this.mapq= rec.getMappingQuality();
		this.samRecord= rec;
		
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
	 * @param samRecord
	 * @param ignoreReadGroup Should the read group be ignored? If true an NA string is returned.
	 * @return
	 */
	private static String getRGtag(SAMRecord samrecord, boolean ignoreReadGroup){
		// Use '\b' as NA since it is not allowed in RG names.
		String rgtag= (samrecord.getAttribute("RG") == null || ignoreReadGroup) ? "\b" : samrecord.getReadGroup().getId();
		return rgtag;
	}

	/* Return object as string suitable to be written to file read
	 * line by line later.
	 * @see java.lang.Object#toString()
	 */
	public String toString(){
		StringBuilder sb= new StringBuilder();
		sb.append(baseQualityScore).append("\t")
		.append(rgtag).append("\t")
		.append(samRecord.getSAMString().trim());
		return sb.toString();
	}
	
	/*  Comparator  */	
	public int compareTo(SAMRecordExt other) {
	    int i = this.samRecord.getReferenceIndex()- other.samRecord.getReferenceIndex();
	    if (i != 0) return i;
	    
	    i = this.samRecord.getUnclippedStart() - other.samRecord.getUnclippedStart();
	    if (i != 0) return i;
	    
	    i = this.samRecord.getUnclippedStart() - other.samRecord.getUnclippedStart();
	    if (i != 0) return i;
	    
	    i= Boolean.valueOf(this.samRecord.getReadNegativeStrandFlag()).compareTo(other.samRecord.getReadNegativeStrandFlag());
	    if (i != 0) return i;
	    
	    i= this.rgtag.compareTo(other.rgtag);
	    if (i != 0) return i;
	    
	    i= other.baseQualityScore - this.baseQualityScore; // DESC: other - this
	    if (i != 0) return i;
	    
	    i= other.samRecord.getMappingQuality() - this.samRecord.getMappingQuality();
	    if (i != 0) return i;
		
	    return i;
	    
	}
	
	/* S E T T E R S  A N D  G E T T E R S */
	
	public SAMRecord getSamRecord() {
		return samRecord;
	}

	public void setSamRecord(SAMRecord samrecord) {
		this.samRecord = samrecord;
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
