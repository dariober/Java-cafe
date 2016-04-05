package bamToBed;

import java.io.File;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class Utils {

	/**
	 * Open sam or bam file given filename or "-" to read from stdin.
	 * @param insam
	 * @param validationStringency See options in ValidationStringency. Typically one of:
	 * 		  ValidationStringency.SILENT, ValidationStringency.LENIENT, ValidationStringency.STRICT 
	 * @return
	 */
	public static SamReader reader(String insam, ValidationStringency validationStringency){
		SamReaderFactory sf = SamReaderFactory.
				makeDefault().
				validationStringency(validationStringency);
		SamReader sam= null;
		if(insam.equals("-")){
			SamInputResource resource= SamInputResource.of(System.in);
			sam= sf.open(resource);			
		} else {
			sam= sf.open(new File(insam));
		}
		return sam;
	}

	/** Convert SAMRecord to bed string */
	public static String SAMRecordToBed(SAMRecord rec) {
		
		// This suffix is to conform to bedtools
		String suffix= "";
		if(rec.getReadPairedFlag()){
			if(rec.getSecondOfPairFlag()){
				suffix= "/2";
			} else if (rec.getFirstOfPairFlag()){
				suffix= "/1";
			}
		}
		StringBuilder sb= new StringBuilder();
		sb.append(rec.getReferenceName());
		sb.append("\t");
		sb.append(rec.getAlignmentStart() - 1);
		sb.append("\t");
		sb.append(rec.getAlignmentEnd());
		sb.append("\t");
		sb.append(rec.getReadName());
		sb.append(suffix);
		sb.append("\t");
		sb.append(rec.getMappingQuality());
		sb.append("\t");
		sb.append(rec.getReadNegativeStrandFlag()? "-" : "+");
		return sb.toString();
	}

	/** Return true if the sam record has to be filtered out. I.e. return true if sam record
	 * DOES NOT have all the bits in requiredFlag or if it has any of the bits in filterFlag 
	 * or the mapq is lower than mapq.*/
	public static boolean filterSamRecord(SAMRecord rec, int requiredFlag, int filterFlag, 
			int mapq) {
		
		if( (rec.getFlags() & requiredFlag) != requiredFlag ){
			return true;
		}
		if( (rec.getFlags() & filterFlag) != 0 ){
			return true;
		} 
		if (rec.getMappingQuality() < mapq){
			return true;
		} 
		return false;
	}
	
}
