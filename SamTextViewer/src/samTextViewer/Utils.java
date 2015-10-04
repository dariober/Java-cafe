package samTextViewer;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.lang3.StringUtils;

import readWriteBAMUtils.ReadWriteBAMUtils;

/**
 * @author berald01
 *
 */
public class Utils {
	/**
	 * Input: 
	  	- Outer list: Lines of output
	  	- Inner list: List of TextRead obj to be printe on same line
	 NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
	 AAAAAAAAAAAA
	  CCCCCCCCCCCC
	              TTTTTTTTTTT
	                           GGGGGGGGGGG
	                                   AAAAAAAA
     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
	 AAAAAAAAAAAA TTTTTTTTTTT  GGGGGGGGGGG       
	  CCCCCCCCCCCC                     AAAAAAAA
	 */
	public static List<List<TextRead>> stackReads(List<TextRead> textReads){
		
		List<List<TextRead>> listOfLines= new ArrayList<List<TextRead>>();
		if(textReads.size() == 0){
			return listOfLines;
		}
		List<TextRead> line= new ArrayList<TextRead>();
		line.add(textReads.get(0)); 
		textReads.remove(0);
		listOfLines.add(line);
		while(true){
			ArrayList<TextRead> trToRemove= new ArrayList<TextRead>();
			// Find a read in input whose start is greater then end of current
			for(int i=0; i < textReads.size(); i++){
				TextRead tr= textReads.get(i);
				if(tr.getTextStart() > line.get(line.size()-1).getTextEnd()+2){ // +2 because we want some space between adjacent reads
					listOfLines.get(listOfLines.size()-1).add(tr); // Append to the last line. 
					trToRemove.add(tr);
					//textReads.remove(i); // If found remove from input 
				}
			} // At the end of the loop you have put in line as many reads as you can. 
			for(TextRead tr : trToRemove){ 
				textReads.remove(textReads.indexOf(tr));
			}
			// Create a new line, add the first textRead in list
			if(textReads.size() > 0){
				line= new ArrayList<TextRead>();
				line.add(textReads.get(0));
				listOfLines.add(line);
				textReads.remove(0);
			} else {
				break;
			}
		}
		return listOfLines;
	}

	/** Prepare a printable string of each output line. 
	 * @param textReads List reads to print out on the same line.
	 * @param fmt Should fomratting be applied? E.g. 2nd in pair bold.
	 * @return
	 */
	public static String linePrinter(List<TextRead> textReads, boolean noFormat){
		StringBuilder sb= new StringBuilder();
		int curPos= 0;
		for(TextRead tr : textReads){
			int nspaces= (tr.getTextStart()-1) - curPos;
			curPos += nspaces;
			sb.append(StringUtils.repeat(" ", nspaces));
				
			for(int i= 0; i < tr.getTextRead().length; i++){
				String fx= (char)tr.getTextRead()[i] + "";
				if(!noFormat){
					// For formatting see http://misc.flogisoft.com/bash/tip_colors_and_formatting
					// and http://stackoverflow.com/questions/5947742/how-to-change-the-output-color-of-echo-in-linux
					if(fx.equals("M") || fx.equals("m")){
						fx= "\033[97;101m" + fx + "\033[0m";
					} else if(fx.equals("U") || fx.equals("u")){
						fx= "\033[97;104m" + fx + "\033[0m";
					}
				}
				sb.append(fx);
				curPos++;
			}
		}
		return sb.toString();
	}
	
	/** Get the coordinates of the first mapped read in bam file and put it
	 * in GenomicCoords obj.
	 * @param bam
	 */
	public static GenomicCoords getStartCoordsOfBAM(String bam){
		SamReader samReader= ReadWriteBAMUtils.reader(bam, ValidationStringency.SILENT);
		String pos= "";
		for(SAMRecord rec: samReader){
			if(!rec.getReadUnmappedFlag()){
				// This is silly: GenomciCoords
				pos= rec.getReferenceName() + ":" + rec.getAlignmentStart() + "-" + rec.getAlignmentStart();  
				//gc.setChrom(rec.getReferenceName());
				//gc.setFrom(rec.getAlignmentStart());
				//gc.setTo(gc.getFrom());
				break;
			}
		}
		GenomicCoords gc= new GenomicCoords(pos, samReader.getFileHeader().getSequenceDictionary());
		try {
			samReader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return gc;
	}
	
	public static GenomicCoords getStartCoordsOfBAM(String bam, String chrom){
		SamReader samReader= ReadWriteBAMUtils.reader(bam, ValidationStringency.LENIENT);
		Iterator<SAMRecord> sam= samReader.query(chrom, 0, 0, false);
		
		GenomicCoords gc= new GenomicCoords(null, null, null, samReader.getFileHeader().getSequenceDictionary());
		while(sam.hasNext()){
			SAMRecord rec= sam.next();
			if(!rec.getReadUnmappedFlag()){
				gc= new GenomicCoords(
					rec.getReferenceName(), 
					rec.getAlignmentStart(),
					rec.getAlignmentStart(),
					samReader.getFileHeader().getSequenceDictionary());
				break;
			}
		}
		try {
			samReader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return gc;		
	}
	
	
	private static int parseStringToIntWithUnits(String x){
		x= x.trim();
		int multiplier= 0;
		if(x.endsWith("k") || x.endsWith("K")){
			multiplier= 1000;
			x= x.substring(0, x.length()-1).trim();
		} else if(x.endsWith("m") || x.endsWith("M")){
			multiplier= 1000000;
			x= x.substring(0, x.length()-1).trim();
		} else if(x.matches("^\\-{0,1}\\d+$") || x.matches("^\\+{0,1}\\d+$")){
			multiplier= 1;
		} else {
			System.err.println("Invalid string to convert to int: " + x);
			System.exit(1);
		}
		int pos= Integer.parseInt(x) * multiplier; 
		return pos;
	}
	
	/**
	 * Parse user to modify the current genomics coordinates in input to new ones
	 * to move.
	 * @param bam 
	 * @return
	 */
	public static String parseConsoleInput(
			String rawInput, GenomicCoords gc){
		
		String region= "";
		String chrom= gc.getChrom();
		Integer from= gc.getFrom();
		Integer to= gc.getTo();
		
		int windowSize= to - from + 1;
		if(rawInput.trim().equals("") || rawInput.trim().equals("f")){				
			from += windowSize; 
			to += windowSize; 
			return chrom + ":" + from + "-" + to;
		} else if(rawInput.trim().equals("b")) {
			from -= windowSize;
			to -= windowSize; 
			return chrom + ":" + from + "-" + to;
		//} else if(rawInput.trim().startsWith("-r")) { // Go to chrom
		//	region= rawInput.trim().substring(2).trim(); // Strip '-r'
		//	return region;
			
		} else if(rawInput.trim().startsWith("+") 
				|| rawInput.trim().startsWith("-") 
				|| Character.isDigit(rawInput.trim().charAt(0))){
			int offset= parseStringToIntWithUnits(rawInput.trim());
			from += offset;
			to += offset;
			return chrom + ":" + from + "-" + to;
		}else if (rawInput.equals("q")) {
			System.exit(0);	
		} else {
			System.err.println("Invalid input for " + rawInput);
		}
		return region;
	}
	
	public static boolean isInteger(String s) {
	 
		s= s.replaceFirst("\\+", ""); // Allow first char to be +
		
		try { 
	        Integer.parseInt(s); 
	    } catch(NumberFormatException e) { 
	        return false; 
	    } catch(NullPointerException e) {
	        return false;
	    }
	    // only got here if we didn't return false
	    return true;
	}
}
