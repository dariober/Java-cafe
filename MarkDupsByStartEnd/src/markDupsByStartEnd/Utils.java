package markDupsByStartEnd;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextTagCodec;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

public class Utils {

	/** Calculates a score for the read which is the sum of scores over Q15. 
	 * Code taken from 
	 * https://github.com/samtools/htsjdk/blob/master/src/java/htsjdk/samtools/DuplicateScoringStrategy.java
	 * */
	public static short getSumOfBaseQualities(final SAMRecord rec) {
		short score = 0;
		for (final byte b : rec.getBaseQualities()) {
			if (b >= 15) score += b;
			}
		return score;
	}

	
	/**
	 * This var is related to samRecordToTabLine().
	 * samRecordToTabLine returns a tab separated string where the last n fields
	 * make up a sam line. OFFSET gives the index of the start of the sam fields.
	 * OFFSET is 0 based, so if the sam record starts at the 8th field, then OFFSET= 7.
	 */
	public static final int OFFSET= 7;

	/**
	 * Turn a SAMRecord to string tab separated useful to define read blocks
	 * @param rec
	 * @param ignoreReadGroup 
	 * @return
	 */
	public static String samRecordToTabLine(SAMRecord rec, boolean ignoreReadGroup){
		/* Create a tab separated file with the following fields.

		MEMO: 1: Order of columns must be consistent with sortTabAndGetOuput()!
		MEMO: 2: Keep OFFSET consistent! 
		*/
		
		String rg= (rec.getAttribute("RG") == null || ignoreReadGroup) ? "\b" : rec.getReadGroup().getId();
		// Integer as= (Integer) ((rec.getAttribute("AS") != null) ? rec.getAttribute("AS") : -1000000); // If AS is NA put a large -ve number.
		// Integer nm= (Integer) rec.getAttribute("NM"); 
		
		StringBuilder sb= new StringBuilder();
		sb.
		append(rec.getContig()).append("\t").         		   // 1.Duplicate defining fields
		append(rec.getUnclippedStart()).append("\t").		   // 2.
		append(rec.getUnclippedEnd()).append("\t").            // 3.
		append(rec.getReadNegativeStrandFlag()).append("\t").  // 4.
		append(rg).append("\t").                               // 5.

		append(getSumOfBaseQualities(rec)).append("\t").       // 6. Discriminatory fields
		append(rec.getMappingQuality()).append("\t").          // 7. 
		
		append(rec.getSAMString());					           // 8. Data (original sam record)
		return sb.toString();
	}
	public static Process sortTabAndGetOuput(String fileToSort) throws IOException, InterruptedException{
		/* Sort tab file by the fields defining the duplicate block,  then by discriminatory fields.
		Now you have reads at the same position, same strand, same RG next to each other,
		in blocks of duplicate reads.
		
		MEMO: Order of columns must be consistent with the one returned by samRecordToTabLine()
		*/
		
		List<String> cmd= new ArrayList<String>();
		cmd.add("sort");
		cmd.add("-t\t");
		cmd.add("-s");
		cmd.add("-S 1G");
		cmd.add("-k1,1");   // chrom
		cmd.add("-k2,2n");  // start
		cmd.add("-k3,3n");  // end
		cmd.add("-k4,4");   // strand
		cmd.add("-k5,5");   // read group
		cmd.add("-k6,6nr"); // sum of base quals
		cmd.add("-k7,7nr"); // mapq
		cmd.add(fileToSort);

		ProcessBuilder pb = new ProcessBuilder(cmd);
		Map<String, String> env = pb.environment();
		env.put("LC_ALL", "C");
		Process p = pb.start();		
		return p;
	}
	
	/**
	 * Convert a string array to a SAMRecord. USeful to make SAMRecords from 
	 * a string which e.g. comes from a tab separated file.
	 * @param array
	 * @param hdr
	 * @return
	 */
	public static SAMRecord arrayToSAMRecord(String[] array, SAMFileHeader hdr){
		SAMRecord rec= new SAMRecord(hdr);
		rec.setReadName(array[0]);
		rec.setFlags(Integer.parseInt(array[1]));
		rec.setReferenceName(array[2]);
		rec.setAlignmentStart(Integer.parseInt(array[3]));
		rec.setMappingQuality(Integer.parseInt(array[4]));
		rec.setCigarString(array[5]);
		rec.setMateReferenceName(array[6]);
		rec.setMateAlignmentStart(Integer.parseInt(array[7]));
		rec.setInferredInsertSize(Integer.parseInt(array[8]));
		rec.setReadString(array[9]);
		rec.setBaseQualityString(array[10]);
		for(int i= 11; i < array.length; i++){
			String tag= array[i];
			Entry<String, Object> tagcode= new TextTagCodec().decode(tag);
			rec.setAttribute(tagcode.getKey(), tagcode.getValue());
		}
		return rec;
	} 
	
	/**
	 * Get a name for a tmp file using "filename" as basename.
	 * Essentially append to filename a suffix making sure that filename+suffix is not
	 * an existing file.
	 * Really you should use File.createTempFile() but I got problems with large tmp files.
	 * @param filename
	 * @return
	 */
	public static String getTmpFilename(String filename){
		String suffix= "markdup.tmp";
		int n= 0;
		String tmpname= filename + "." + n + "." + suffix;
		while( new File(tmpname).isFile() ){
			n++;
			tmpname= filename + "." + n + "." + suffix;
		}
		return tmpname;
	}
		
}
