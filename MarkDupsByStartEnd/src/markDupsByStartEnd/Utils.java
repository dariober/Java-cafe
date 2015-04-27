package markDupsByStartEnd;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextTagCodec;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

public class Utils {

	/**
	 * Turn a SAMRecord to string tab separated useful to define read blocks
	 * @param rec
	 * @return
	 */
	public static String samRecordToTabLine(SAMRecord rec){
		/* Create a tab separated file with the following fields.

		MEMO: Order of columns must be consistent with sortTabAndGetOuput()
		*/

		String rg= (rec.getAttribute("RG") != null) ? rec.getReadGroup().getId() : "";
		Integer as= (Integer) ((rec.getAttribute("AS") != null) ? rec.getAttribute("AS") : -1000000);
		Integer nm= (Integer) rec.getAttribute("NM"); 
		
		StringBuilder sb= new StringBuilder();
		sb.
		append(rec.getContig()).append("\t").         // Duplicate defining fields
		append(rec.getUnclippedStart()).append("\t").
		append(rec.getUnclippedEnd()).append("\t").
		append(rec.getReadNegativeStrandFlag()).append("\t").
		append(rg).append("\t").
		
		append(rec.getMappingQuality()).append("\t"). // Discriminatory fields
		append(as).append("\t").
		append(nm).append("\t").
		
		append(rec.getSAMString());					  // Data (original sam record)
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
		cmd.add("-s");
		cmd.add("-S 1G");
		cmd.add("-k1,1");   // chrom
		cmd.add("-k2,2n");  // start
		cmd.add("-k3,3n");  // end
		cmd.add("-k4,4");   // strand
		cmd.add("-k5,5");   // read group
		cmd.add("-k6,6nr"); // mapq
		cmd.add("-k7,7nr"); // Alignment score AS
		cmd.add("-k8,8n");  // Edit distance NM
		cmd.add(fileToSort);
		
		// ** Uncomment for debugging **
		//StringBuilder sb= new StringBuilder();
		//for(String x : cmd){
		//	sb.append(x + " ");
		//}
		// System.err.println(sb.toString());
		ProcessBuilder pb = new ProcessBuilder(cmd);
		Map<String, String> env = pb.environment();
		env.put("LC_ALL", "C");
		Process p = pb.start();
		return p;
	}

	// UNUSED
	public static void checkExitSortTabAndGetOuput(Process p) throws IOException{
		
		if(p.exitValue() != 0){
			System.err.println("Exited with code: " + p.exitValue());
			BufferedReader er = new BufferedReader(
		            new InputStreamReader(p.getErrorStream()));
			String line = null;
			while ((line = er.readLine()) != null) {
				System.err.println(line);
			}
			System.exit(1);
		} 
		
	}
	
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
		
}
