package samTextViewer;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Set;

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
					String fmt= "\033["; // Start format 
					if(tr.rec.getReadPairedFlag() && tr.rec.getSecondOfPairFlag()){
						fmt += "4;"; // Underline 2nd n pair
					}					
					if(fx.toUpperCase().equals("M")){
						fmt += "97;101";
					} else if(fx.toUpperCase().equals("U")){
						fmt += "97;104";
					} else if(fx.toUpperCase().equals("A")){
						fmt += "1;107;34";  //1: bold; 107: white bg
					} else if(fx.toUpperCase().equals("C")) {
						fmt += "1;107;31";
					} else if(fx.toUpperCase().equals("G")) {
						fmt += "1;107;32";
					} else if(fx.toUpperCase().equals("T")) {
						fmt += "1;107;33";
					}
					// The formatted string will look like `echo -e "\033[4;1;107;31mACTGnnnnnACTG\033[0m"`
					fx= fmt + "m" + fx + "\033[0m"; // Clear all formatting
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
		
		GenomicCoords gc= new GenomicCoords(null, samReader.getFileHeader().getSequenceDictionary());
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
		if(rawInput.trim().equals("ff")){				
			from += windowSize; 
			to += windowSize; 
			return chrom + ":" + from + "-" + to;
		} else if(rawInput.trim().equals("bb")) {
			from -= windowSize;
			to -= windowSize; 
			return chrom + ":" + from + "-" + to;
		} else if(rawInput.trim().equals("f")){
			int step= (int)Math.round(windowSize / 10d);
			step= (step == 0) ? 1 : step;
			from += step; 
			to += step;
			return chrom + ":" + from + "-" + to;
		} else if(rawInput.trim().equals("b")){
			int step= (int)Math.round(windowSize / 10d);
			step= (step == 0) ? 1 : step;
			from -= step; 
			to -= step;
			return chrom + ":" + from + "-" + to;
		} else if(rawInput.trim().startsWith(":")) { // You might want to be more specific than just startsWith(:)
			String pos= rawInput.trim().replaceFirst(":", "");
			Integer.parseInt(pos); // Check you actually got an int.
			return chrom + ":" + pos;
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

	/**
	 * Compress the input list of ints in a list of length nwinds by 
	 * grouping elements and taking the mean of each group (or other summary stat).
	 * 
	 * **NB**: If you change the splitting algorithm here check also Ruler_TO_BE_DEPRECTED class is
	 * consistent.
	 * 
	 * @param ints
	 * @param nwinds
	 * @return Map where key is index (0-based) of first element in original input 
	 * and value is averaged of grouped elements. The keys can be used as ruler for the coverage track.  
	 */
	/*public static LinkedHashMap<Integer, Double> compressNumericListTO_BE_DEPRECATED(List<Double> ints, int nwinds){
		
		int grpSize= (int) Math.round(((double)ints.size() / nwinds));
		// After round() you get a remainder which goes in the last bin
		LinkedHashMap<Integer, Double> zlist= new LinkedHashMap<Integer, Double>();
		List<Double> sublist= new ArrayList<Double>();
		// int i= 0;
		int at= 0;
		for(int i= 0; i < ints.size(); i++){
			sublist.add(ints.get(i));

			if(sublist.size() == grpSize || grpSize < 1){ // < 1 is for num. of windows >num. elements. 
														  // So no compression done 
				double avg= calculateAverage(sublist);
				zlist.put(at, avg);
				sublist.clear();
				at= i+1;
			
			}
		}
		if(sublist.size() > 0){
			double avg= Math.round(calculateAverage(sublist));
			zlist.put(at, avg);			
		}
		return zlist;
	} */

	
	/**
	 * Average of ints in array x. Adapted from:
	 * http://stackoverflow.com/questions/10791568/calculating-average-of-an-array-list
	 * null values are ignored, like R mean(..., na.rm= TRUE). 
	 * Returns Float.NaN if input list is empty or only nulls. You can check for Float.NaN
	 * with Float.isNaN(x); 
	 * @param marks
	 * @return
	 */
	public static Double calculateAverage(List<Double> list) {
		double sum = 0;
		long  N= 0;  
		if(!list.isEmpty()) {
			for (Double z : list) {
				if(z  != null && !Double.isNaN(z)){
					sum += z;
					N++;
				}
			}
			return (double)sum / N;
		}
		return Double.NaN;
	}

	/**
	 * Naive search to get the index position of the value in list closest to a given value.
	 * 
	 * @param genomePos
	 * @param mapping
	 * @return
	 */
	public static int getIndexOfclosestValue(double genomePos, List<Double> mapping){
		double bestDiff= Integer.MAX_VALUE;
		int closest= -1;
		for(int idx= 0; idx < mapping.size(); idx++){ 
			// Iterate through entire list to find closest position on screen, it's a bit wasteful since
			// the list is ordered, but it's ok.
			double candidate= mapping.get(idx);
			double diff= Math.abs(genomePos - candidate);
			if(diff < bestDiff){
				closest= idx;
				bestDiff= diff;
			}
		}
		if(closest < 0){
			System.err.println("Invalid index position.");
			System.exit(1);
		}
		return closest;
	}
		
}
