package samTextViewer;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.lang3.StringUtils;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

/** Class holding features of a read necessary to represent the read in text
 * format.
 * Missing/Todo: Insertions to the reference are not visible.
 * @author berald01
          10        20        30        40        
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
AACCTTGGCC---------------------------------------
-AACCTTGGCC--------------------------------------
--AACCTTGGCC-------------------------------------
-------------AACCTTGGCC--------------------------
--------------AA----CCTT-------------------------
 */
public class TextRead {
	
	/** Char to represent deletions from the reference. I.e. gaps in the read */
	final private byte DEL= '-';
	/** Char to represent region skip. I.e. gaps in the read */
	final private byte N= '_';
	
	/** Start position of the read in window coordinates. 1-based. If in genomic
	 * coords read starts aligning at pos 100 and window span is 100:150, then textStart= 1.*/
	private int textStart;

	/** End position of the read in window coordinates. 1-based. If in genomic
	 * coords read ends pos 150 and window span is 100:150, then textEnd= 1. */
	private int textEnd;
	
	// The read represented as text: The actual bases as found in SAMRecord, 
	// and the bases converted to '.' and ',' as appropriate.
	private byte[] textRead;  // Actual read bases or DEL if gaps 
	public SAMRecord rec;
	
	/*    C o n s t r u c t o r s    */
		
	public TextRead(SAMRecord rec, TextWindow textWindow){
		// At least some of the read must be in the window
		if(rec.getAlignmentStart() > textWindow.getTo()){
			System.err.println("Alignment starts beyond text window!");
			System.err.println(rec.getSAMString());
			System.err.println("Aln starts: " +  rec.getAlignmentStart());
			System.err.println("Aln ends: " +  rec.getAlignmentEnd());
			System.err.println("Window: " + textWindow.toString());
			System.exit(1);
		}
		if(rec.getAlignmentEnd() < textWindow.getFrom()){
			System.err.println("Alignment ends before text window!");
			System.err.println(rec.getSAMString());
			System.err.println("Aln starts: " +  rec.getAlignmentStart());
			System.err.println("Aln ends: " +  rec.getAlignmentEnd());
			System.err.println("Window: " + textWindow.toString());
			System.exit(1);
		}		
		this.textStart= getTextStart(rec, textWindow);
		this.textEnd= getTextEnd(rec, textWindow);
		this.textRead= getTextRead(rec, textWindow);
		this.rec= rec;
	}
	/*       M e t h o d s       */
	
	public int getTextStart(SAMRecord rec, TextWindow textWindow){
		textStart= rec.getAlignmentStart() - textWindow.getFrom() + 1;
		if(textStart <= 0){
			textStart= 1;
		}
		return textStart;
	}
	
	public int getTextEnd(SAMRecord rec, TextWindow textWindow){
		int windowSize= textWindow.getTo() - textWindow.getFrom();
		if(rec.getAlignmentEnd() > textWindow.getTo()){
			textEnd= windowSize;
		} else {
			textEnd= rec.getAlignmentEnd() - textWindow.getFrom();
		}
		return textEnd;
	}
	/** Get a representation of the read as it appears aligned to the reference. 
	 * I.e. clipped ends omitted and deletions appearing as gaps (empty byte).
	 * Only the portion contained between the genomic coords from:to is returned.
	 * @return
	 */
	private byte[] getTextRead(SAMRecord rec, TextWindow textWindow) {

		// Accumulate here the read bases inside the window 
		ArrayList<Byte> textRead= new ArrayList<Byte>();
		byte[] readBases= rec.getReadBases();
		// Walk along the aligned read and append bases to textRead as long as
		// the genomic position of the base is inside the genomic coords of the window
		int curBaseGenomicPos= rec.getAlignmentStart();
		int curBaseReadPos= 0; // Position on read. Start from zero walk along the read
		List<CigarElement> cigarEls= rec.getCigar().getCigarElements();
		for(CigarElement el : cigarEls){
			if(el.getOperator() == CigarOperator.M || el.getOperator() == CigarOperator.EQ || el.getOperator() == CigarOperator.X){
				for(int i= 0; i < el.getLength(); i++){
					if(curBaseGenomicPos >= textWindow.getFrom() && curBaseGenomicPos <= textWindow.getTo()){
						// If base is inside window:
						if(readBases.length > 0){
							textRead.add(readBases[curBaseReadPos]);
						} else { // If sam record has no read seq stored put N
							textRead.add((byte)'N');
						}
					}
					curBaseGenomicPos++; // M consumes read and ref bases. So increment them
					curBaseReadPos++;
				}
			} else if(el.getOperator() == CigarOperator.D || el.getOperator() == CigarOperator.N){
				for(int i= 0; i < el.getLength(); i++){
					if(curBaseGenomicPos >= textWindow.getFrom() && curBaseGenomicPos <= textWindow.getTo()){
						if(el.getOperator() == CigarOperator.D){
							textRead.add(DEL);
						} else if(el.getOperator() == CigarOperator.N){ 
							textRead.add(N);
						} else {
							System.err.println("Unexpected operator"); System.exit(1);
						}
					}
					curBaseGenomicPos++;
				}
			} else if(el.getOperator() == CigarOperator.I) {
				curBaseReadPos += el.getLength(); // Insertions in the reference are missed
			} else if(el.getOperator() == CigarOperator.S){
				curBaseReadPos += el.getLength();
			} else if(el.getOperator() == CigarOperator.H){
				// Nothing to do
			} else if(el.getOperator() == CigarOperator.P){
				// Nothing to do: NOT SURE is is correct to just ignore padding!
			} else {
				System.err.println("Unexpected operator in cigar string for record\n" + rec.getSAMString()); System.exit(1);
			}
		}
		byte[] textReadArr= new byte[textRead.size()];
		for(int i=0; i < textRead.size(); i++){
			textReadArr[i]= textRead.get(i);
		}
		return textReadArr;
	}
	
	/** Convert textRead, the actual bases found in sam, to represent match, mismatch
	 * and strandness.
	 * @param refSeq The reference sequence spanning the window.
	 * @return
	 */
	public byte[] getTextChars(byte[] refSeq) {
				
		byte[] readChars= new byte[this.textRead.length];
		int posOnRead= 0;
		for(int i= this.textStart - 1; i <= this.textEnd; i++){
			byte base= (byte) Character.toUpperCase( this.textRead[posOnRead] );
			byte ref= (byte) Character.toUpperCase(refSeq[i]);
			if( base == ref){
				if(this.rec.getReadNegativeStrandFlag()){
					readChars[posOnRead++]= ',';
				} else {
					readChars[posOnRead++]= '.';
				}
			} else {
				if(this.rec.getReadNegativeStrandFlag()){
					readChars[posOnRead++]= (byte)Character.toLowerCase(base);
				} else {
					readChars[posOnRead++]= base;
				}				
			}
		}
		return readChars;
	}

	/** Memo: You compare the reads with these bases on the reference:
	   1st |  2nd
	-------+------
	+ve  C |  G
	-------+------
 	-ve  g |  c
 	
	 * @param refSeq
	 * @return
	 */
	public byte[] convertTextToBS(byte[] refSeq){
				
		byte[] readChars= this.getTextChars(refSeq);
		byte[] textBS= new byte[readChars.length];
		
		for(int i=0; i < refSeq.length; i++){
			// Now the reference is completely uppercase
			refSeq[i]= (byte) Character.toUpperCase(refSeq[i]);
		}
		
		// For convenience extract flags from sam record
		boolean isSecondOfPair= false;
		if(rec.getReadPairedFlag() && rec.getSecondOfPairFlag()){
			isSecondOfPair= true;
		}
		boolean isForwardStrand= !this.rec.getReadNegativeStrandFlag();
		
		for(int i= 0; i < textBS.length; i++){
			byte ref= refSeq[i + this.getTextStart() - 1];
			byte read= readChars[i];
			if( ( isForwardStrand && !isSecondOfPair ) || ( !isForwardStrand && isSecondOfPair )){
				// Look for C on the reference
				if(isForwardStrand){ // +strand, first in pair or unpaired
					if(ref == 'C' && read == '.'){
						textBS[i]= 'M';
					} else if(ref == 'C' && read == 'T'){
						textBS[i]= 'U';
					} else {
						textBS[i]= read;
					}
				} else { // -ve strand, 2nd in pair
					// Look for c=',' -> m; 't' -> u
					if(ref == 'C' && read == ','){
						textBS[i]= 'm';
					} else if(ref == 'C' && read == 't'){
						textBS[i]= 'u';
					} else {
						textBS[i]= read;
					}
				}
			} else if( ( !isForwardStrand && !isSecondOfPair ) || ( isForwardStrand && isSecondOfPair )){
				// Look for G on the reference
				if(!isForwardStrand){ // -ve strand; first in pair or unpaired
					if(ref == 'G' && read == ','){
						textBS[i]= 'm';
					} else if(ref == 'G' && read == 'a'){
						textBS[i]= 'u';
					} else {
						textBS[i]= read;
					}
				} else { // -ve strand, 2nd in pair
					if(ref == 'G' && read == '.'){
						textBS[i]= 'M';
					} else if(ref == 'G' && read == 'A'){
						textBS[i]= 'U';
					} else {
						textBS[i]= read;
					}
				}
			}
		}
		return textBS;
	}
	
	/** Combine the text representation of the read with the starting pos on the
	 * window to produce a printable string of the read.
	 * @return
	 */
	public String getPrintableText(){
		StringBuilder sb= new StringBuilder();
		if(this.textStart > 1){
			sb.append(StringUtils.repeat(" ", this.textStart - 1));
		}
		sb.append(new String(this.textRead));
		return sb.toString();
	}

	public String getPrintableText(byte[] refSeq){
		StringBuilder sb= new StringBuilder();
		byte[] textChars= this.getTextChars(refSeq);
		if(this.textStart > 1){
			sb.append(StringUtils.repeat(" ", this.textStart - 1));
		}
		sb.append(new String(textChars));
		return sb.toString();		
	}
	
	public String toString(){
		StringBuilder sb= new StringBuilder();
		sb.append("Text start: " + this.textStart + "\n");
		sb.append("Text end: " + this.textEnd + "\n");
		sb.append("Text read: [" + new String(this.textRead) + "]\n");
		sb.append("Printable read: [" + new String(this.getPrintableText()) + "]");
		return sb.toString();
	}
	
	/*      S e t t e r s   and   G e t t e r s     */
	public int getTextStart() {
		return textStart;
	}

	public void setTextStart(int textStart) {
		this.textStart = textStart;
	}

	public int getTextEnd() {
		return textEnd;
	}

	public void setTextEnd(int textEnd) {
		this.textEnd = textEnd;
	}

	public byte[] getTextRead() {
		return textRead;
	}

	public void setTextRead(byte[] textRead) {
		this.textRead = textRead;
	}

}
