package vmatchToSam;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import exceptions.InvalidVmatchRecordException;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.ValidationStringency;
import utils.Utils;

public class VmatchRecord {

	public static String DEFAULT_READ_GROUP= "NA";
	
	private int referenceAlignmentLength;
	private String referenceNameOrIndex;
	private int referenceAlignmentStart; // 1-based
	private Object strand;
	private int queryAlignmentLength; // N of bases of the query aligned to ref. Gaps don't count.
	private String queryNameOrIndex;
	private int queryAlignmentStart;
	private int editDistance;
	private float evalue;
	private int alignmentScore;
	private float pctIdentity;
	private StringBuilder recordString= new StringBuilder();
	private String alignedReferenceSequence;
	private String alignedQuerySequence;
	
	/*      C O N S T R U C T O R S      */

	public VmatchRecord() {
		
	}
	
	/**
	 * Read the next vmatch record from BufferedReader. 
	 * @throws IOException 
	 * @throws InvalidVmatchRecordException 
	 * */
	public VmatchRecord(BufferedReader br) throws IOException, InvalidVmatchRecordException {
		int nConsecutiveBreaks= 0;
		boolean alnStatsNeeded= true;
		StringBuilder referenceSequence= new StringBuilder();
		StringBuilder querySequence= new StringBuilder();
		while(nConsecutiveBreaks < 2){
			String line= br.readLine().trim();
			
			if(line.startsWith("#")){
				continue;
			}
			if(line.isEmpty()){
				nConsecutiveBreaks += 1;
			} else {
				nConsecutiveBreaks= 0;
			}
			this.recordString.append(line + "\n");
			if(alnStatsNeeded){
				this.parseAlnStatsString(line);
				alnStatsNeeded= false;
			} else if(line.startsWith("Sbjct:") && referenceSequence.length() == querySequence.length()) {
				// Alternate between appending sequence to reference and query
				referenceSequence.append(line.replaceFirst(".+? ", "").replaceAll(" .*", ""));
			} else if((line.startsWith("Sbjct:") || line.startsWith("Query:")) && querySequence.length() < referenceSequence.length()){
				querySequence.append(line.replaceFirst(".+? ", "").replaceAll(" .*", ""));
			} else {
				continue;
			}
		}
		
		this.alignedReferenceSequence= referenceSequence.toString();
		this.alignedQuerySequence= querySequence.toString();

		// Some checks
		// ===========
		if(this.alignedReferenceSequence.length() != this.alignedQuerySequence.length()){
			System.err.println("Reference: " + this.alignedReferenceSequence);
			System.err.println("Query:     " + this.alignedQuerySequence);
			throw new InvalidVmatchRecordException();
		}
		if(this.alignedReferenceSequence.replaceAll("-", "").length() != this.referenceAlignmentLength){
			System.err.println("Ref seq length: " + this.alignedReferenceSequence.length());
			System.err.println("Ref aln len:    " + this.referenceAlignmentLength);
			throw new InvalidVmatchRecordException();
		}
		if(this.alignedQuerySequence.replaceAll("-", "").length() != this.queryAlignmentLength){
			System.err.println("Query seq length: " + this.alignedQuerySequence.length());
			System.err.println("Query aln len:    " + this.queryAlignmentLength);
			throw new InvalidVmatchRecordException();
		}
	}

	/* M E T H O D S */
	
	/** Populate vmatch record by parsing the line containing the 
	 * alignment stats. Typically this looks like:
	 * 25   ref_1   2   D 27   r1   0   2    7.44e-10  46    92.5
	 * @throws InvalidVmatchRecordException 
	 * */
	private void parseAlnStatsString(String alnStats) throws InvalidVmatchRecordException {
		String[] stats= alnStats.trim().split(" +");
	
		this.referenceAlignmentLength= Integer.parseInt(stats[0]);
		this.referenceNameOrIndex= stats[1];
		this.referenceAlignmentStart= Integer.parseInt(stats[2]);
		this.strand= stats[3].equals("D") ? "+" : stats[3].equals("P") ? "-" : null;
		if(this.strand == null){
			throw new InvalidVmatchRecordException(); 
		}
		this.queryAlignmentLength= Integer.parseInt(stats[4]);
		this.queryNameOrIndex= stats[5];
		this.queryAlignmentStart= Integer.parseInt(stats[6]);
		this.editDistance= Integer.parseInt(stats[7]);
		this.evalue= Float.parseFloat(stats[8]);
		this.alignmentScore= Integer.parseInt(stats[9]); 
		this.pctIdentity= Float.parseFloat(stats[10]);
	}

	public SAMRecord getSAMRecord() {
		SAMRecord samRec= new SAMRecord(null);
		samRec.setReadName(this.queryNameOrIndex);
		if(this.strand.equals("-")){
			samRec.setReadNegativeStrandFlag(true);
		}
		// Applying refSeqNameForSam() shouldn't be necessary as vmatch already strips the blank and what follows  
		samRec.setReferenceName(Utils.refSeqNameForSam(this.referenceNameOrIndex));
		samRec.setAlignmentStart(this.referenceAlignmentStart + 1);
		samRec.setMappingQuality(255);
		samRec.setCigar(this.makeCigarFromAlnStrings());
		samRec.setReadBases(this.alignedQuerySequence.replaceAll("-", "").getBytes());

		// Attributes
		samRec.setAttribute("NM", this.editDistance);
		samRec.setAttribute("AS", this.alignmentScore);
		samRec.setAttribute("XE", this.evalue);
		samRec.setAttribute("XP", this.pctIdentity);
		samRec.setAttribute("RG", DEFAULT_READ_GROUP);
		// I don't think the cigar is validated against the sequence though...  
		samRec.setValidationStringency(ValidationStringency.STRICT);
		samRec.validateCigar(-1);
		return samRec;
	}

	/**
	 * Compare aligned reference and query strings to get cigar
	 *  Memo:
	 *  --ACTG--ACTG---AATTTTAAA : Ref
	 *  TTACTGAAACTGCCCAA---T--- : Qry
	 *  SSMMMMIIMMMMIIIMMDDDMHHH
	 * */
	private Cigar makeCigarFromAlnStrings() {

		// Count number of leading '-' in aligned query: This is the start H operator
		// Count number of leading '-' in aligned reference: This is the start S operator
		int startH= this.alignedQuerySequence.length() - this.alignedQuerySequence.replaceAll("^-+", "").length();
		int startS= this.alignedReferenceSequence.length() - this.alignedReferenceSequence.replaceAll("^-+", "").length();
		
		// Count number of trailing '-' in aligned query: This is the end H operator
		// Count number of trailing '-' in aligned reference: This is the end S operator
		int endH= this.alignedQuerySequence.length() - this.alignedQuerySequence.replaceAll("-+$", "").length();
		int endS= this.alignedReferenceSequence.length() - this.alignedReferenceSequence.replaceAll("-+$", "").length();
		
		// Start and end indexes of the matching part of the alignment
		int matchStart= Math.max(startH, startS);
		int matchEnd= this.alignedQuerySequence.length() - Math.max(endH, endS);
			
		List<CigarOperator> cigarOps= new ArrayList<CigarOperator>();
		List<Integer> cigarCnt= new ArrayList<Integer>();
		// Start cigar lists with soft and hard clips
		cigarOps.add(0, CigarOperator.HARD_CLIP);
		cigarCnt.add(0, this.queryAlignmentStart + startH);
		cigarOps.add(1, CigarOperator.SOFT_CLIP);
		cigarCnt.add(1, startS);
		for(int i= matchStart; i < matchEnd; i++){
			char ref= this.alignedReferenceSequence.charAt(i);
			char qry= this.alignedQuerySequence.charAt(i);
			CigarOperator op;
			if(ref == '-'){
				op= CigarOperator.INSERTION;
			} else if(qry == '-'){
				op= CigarOperator.DELETION;
			} else if(ref == qry){
				op= CigarOperator.EQ;
			} else {
				op= CigarOperator.X;
			}
			// Decide whether the operator increments the previous one or it starts a new one
			if(cigarOps.get(cigarOps.size()-1).equals(op)){
				int cur= cigarCnt.get(cigarCnt.size()-1);
				cigarCnt.set(cigarCnt.size()-1, cur+1);
			} else {
				cigarOps.add(op);
				cigarCnt.add(1);
			}
		}
		// Append to cigar lists the soft
		// End hard clip should be the query length - N. aligned bases. However, the query length is not 
		// given in vmatch output. 
		cigarOps.add(CigarOperator.HARD_CLIP);
		cigarCnt.add(this.queryAlignmentStart + endH);
		cigarOps.add(CigarOperator.SOFT_CLIP);
		cigarCnt.add(endS);

		// Build cigar
		Cigar cigar= new Cigar();
		for(int i= 0; i < cigarOps.size(); i++){
			if(cigarCnt.get(i) > 0){
				CigarElement el= new CigarElement(cigarCnt.get(i), cigarOps.get(i));
				cigar.add(el);
			}
		}
		return cigar;
	}

	/** Return the lines read from BufferedReader */
	public String toString(){
		return this.recordString.toString().trim();
	}	
	
	/*   S E T T E R S   AND   G E T T E R S   */
	
	public int getReferenceAlignmentLength() {
		return referenceAlignmentLength;
	}

	public String getReferenceNameOrIndex() {
		return referenceNameOrIndex;
	}

	public int getReferenceAlignmentStart() {
		return referenceAlignmentStart;
	}

	public Object getStrand() {
		return strand;
	}

	public int getQueryAlignmentLength() {
		return queryAlignmentLength;
	}

	public String getQueryNameOrIndex() {
		return queryNameOrIndex;
	}

	public int getQueryAlignmentStart() {
		return queryAlignmentStart;
	}

	public int getEditDistance() {
		return editDistance;
	}

	public double getEvalue() {
		return evalue;
	}

	public double getAlignmentScore() {
		return alignmentScore;
	}

	public double getPctIdentity() {
		return pctIdentity;
	}
	
	public String getAlignedReferenceSequence() {
		return alignedReferenceSequence;
	}
	
	public String getAlignedQuerySequence() {
		return alignedQuerySequence;
	}

}
