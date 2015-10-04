package samTextViewer;

import htsjdk.samtools.SAMRecord;

public class SamTextViewer {

	private String readWithSpacing;
	private String readAligned;
	
	final public char DEL= '-';
	
	/*       C O N S T R UC T O R S       */
	public SamTextViewer(){
		
	};
	
	public SamTextViewer(SAMRecord rec, String ref, int from, int to){
		String readWithSpacing= readWithSpacing(rec, from, to).toUpperCase();
		ref= ref.toUpperCase();
		assert readWithSpacing.length() == ref.length();

		// Convert read to read aligned to reference with appropriate encoding
		char[] readAlignedArr= new char[readWithSpacing.length()];
		for(int i= 0; i < readWithSpacing.length(); i++){
			char readChar= readWithSpacing.charAt(i);
			char refChar= ref.charAt(i);
			if(readChar == ' '){
				readAlignedArr[i]= '-';
			} else if(readChar == refChar){
				readAlignedArr[i]= '.';
			} else {
				readAlignedArr[i]= readChar;
			}
		}
		String readAligned= new String(readAlignedArr);
		if(rec.getReadNegativeStrandFlag()){
			readAligned= readAligned.toLowerCase().replace('.', ',');
		}
		this.setReadAligned(readAligned);
		this.readWithSpacing= readWithSpacing;
	};
	
	private static String readWithSpacing(SAMRecord rec, int from, int to) {

		char[] spacing= new char[(to - from)]; 
		for(int i= 0; i < spacing.length; i++){
			spacing[i]= ' ';
		}
		int offset= rec.getAlignmentStart() - from;

		// MEMO: samtools is 1-based. First base of the read has index 1.
		int readAlnStart= (rec.getAlignmentStart() - rec.getUnclippedStart()) + 1;

		// Number of read bases aligned: Read length minus left and right soft clipped bases.
		int readAlignedLen= rec.getReadLength() - (readAlnStart - 1) - (rec.getUnclippedEnd() - rec.getAlignmentEnd());
				
		for(int i= 0; i < readAlignedLen; i++){
			int readOffset= i + readAlnStart; // Skip left soft clipped bases
			int refPos= rec.getReferencePositionAtReadPosition(readOffset);
			char readBase= (char)rec.getReadBases()[readOffset-1];
			int spaceOffset= offset + (refPos - rec.getAlignmentStart());
			if(spaceOffset >= 0 && spaceOffset < spacing.length){
				spacing[spaceOffset]= readBase;
			}
		}
		return new String(spacing);
	}

	/**
	 * @return
	 */
	//private String fillDeletions(String ){
	//	
	//}
	
	/*    S E T T E R S   AND   G E T T E R S   */
	
	public String getReadWithSpacing() {
		return readWithSpacing;
	}

	public void setReadWithSpacing(String readWithSpacing) {
		this.readWithSpacing = readWithSpacing;
	}

	public String getReadAligned() {
		return readAligned;
	}

	public void setReadAligned(String readAligned) {
		this.readAligned = readAligned;
	}

}
