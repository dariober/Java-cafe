package samTextViewer;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.tribble.readers.TabixReader;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.zip.GZIPInputStream;

import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.tdf.TDFReader;

import tracks.IntervalFeatureSet;
import tracks.IntervalFeature;
import tracks.TrackFormat;

/**
 * @author berald01
 *
 */
public class Utils {
	
	/** Return true if fileName has a valid tabix index. 
	 * @throws IOException 
	 * */
	public static boolean hasTabixIndex(String fileName) throws IOException{
		try{
			TabixReader tabixReader= new TabixReader(fileName);
			tabixReader.readLine();
			tabixReader.close();
			return true;
		} catch (Exception e){
			return false;
		}
	}
	
	/** Get the first chrom string from first line of input file. As you add support for more filetypes you should update 
	 * this function. This method is very dirty and shouldn't be trusted 100% */
	public static String initRegionFromFile(String x) throws IOException{
		String region= "";
		TrackFormat fmt= Utils.getFileTypeFromName(x); 
		if(fmt.equals(TrackFormat.BAM)){
			SamReaderFactory srf=SamReaderFactory.make();
			srf.validationStringency(ValidationStringency.SILENT);
			SamReader samReader = srf.open(new File(x));
			region= samReader.getFileHeader().getSequence(0).getSequenceName();
			samReader.close();
			return region;
		} else if(fmt.equals(TrackFormat.BIGWIG)){
			BBFileReader reader= new BBFileReader(x);
			region= reader.getChromosomeNames().get(0);
			reader.close();
			return region;
		} else if(fmt.equals(TrackFormat.TDF)){
			Iterator<String> iter = TDFReader.getReader(x).getChromosomeNames().iterator();
			while(iter.hasNext()){
				region= iter.next();
				if(!region.equals("All")){
					return region;
				}
			} 
			System.err.println("Cannot initialize from " + x);
			throw new RuntimeException();
		} else {
			// Input file appears to be a generic interval file. We expect chrom to be in column 1
			BufferedReader br;
			if(x.toLowerCase().endsWith(".gz")){
				InputStream fileStream = new FileInputStream(x);
				GZIPInputStream gzipStream = new GZIPInputStream(fileStream);
				Reader decoder = new InputStreamReader(gzipStream, "UTF-8");
				br = new BufferedReader(decoder);
			} else {
				br = new BufferedReader(new FileReader(x));
			}
			String line;
			while ((line = br.readLine()) != null){
				line= line.trim();
				if(line.startsWith("#") || line.isEmpty()){
					continue;
				}
				IntervalFeature feature= new IntervalFeature(line, fmt);
				region= feature.getChrom() + ":" + feature.getFrom(); 
				br.close();
				return region;
			}
		} 
		System.err.println("Cannot initialize from " + x);
		throw new RuntimeException();
	}
	
	public static boolean bamHasIndex(String bam) throws IOException{

		SamReaderFactory srf=SamReaderFactory.make();
		srf.validationStringency(ValidationStringency.SILENT);
		SamReader samReader = srf.open(new File(bam));
		boolean hasIndex= samReader.hasIndex();
		samReader.close();
		return hasIndex;
		
	}
	
	/*public static GenomicCoords findNextStringOnFile(String string, String filename, GenomicCoords curGc,
			Map<String, IntervalFeatureSet> intervalFiles ) throws InvalidGenomicCoordsException, IOException{

		String chosenFn= "";
		if(filename.isEmpty() && intervalFiles.size() == 1){ // Only one file to chose from: Get that one
			chosenFn= new ArrayList<String>(intervalFiles.keySet()).get(0);
		} else {
			// Try to match file perfectly as it was added from cli, including path if any
			for(String fn : intervalFiles.keySet()){
				if(fn.equals(filename)){
					chosenFn = fn;
				}
			}
			if(chosenFn.isEmpty()){
				// Or try to match only file name.
				for(String fn : intervalFiles.keySet()){ // Do not look for a perfect match since the original input might contain path. 
					String onlyName= new File(fn).getName();
					if(onlyName.equals(filename)){
						chosenFn = fn;
					}
				}
			}
		}
		if(chosenFn.isEmpty()){
			System.err.println("File " + filename + " not found in file set\n" + intervalFiles.keySet());
			return curGc;
		}
		return intervalFiles.get(chosenFn).findNextString(curGc, string);
	} */
	
	public static TrackFormat getFileTypeFromName(String fileName){
		fileName= fileName.toLowerCase();
		
		if(    fileName.endsWith(".bed") 
		    || fileName.endsWith(".bed.gz") 
		    || fileName.endsWith(".bed.gz.tbi")){
			return TrackFormat.BED;
		} else if( fileName.endsWith(".gtf") 
				|| fileName.endsWith(".gtf.gz")
				|| fileName.endsWith(".gtf.gz.tbi")
				|| fileName.endsWith(".gff") 
				|| fileName.endsWith(".gff.gz") 
				|| fileName.endsWith(".gff.gz.tbi")
				|| fileName.endsWith(".gff3")
				|| fileName.endsWith(".gff3.gz") 
				|| fileName.endsWith(".gff3.gz.tbi")){
			return TrackFormat.GFF;
		} else if(fileName.endsWith(".bam") || fileName.endsWith(".cram")){
			return TrackFormat.BAM;
		} else if(fileName.endsWith(".bigwig") || fileName.endsWith(".bw")) {
			return TrackFormat.BIGWIG;
		} else if(fileName.endsWith(".tdf")) {
			return TrackFormat.TDF;
		} else if(fileName.endsWith(".bedgraph.gz") || fileName.endsWith(".bedgraph")) {
			return TrackFormat.BEDGRAPH;
		} else {
			// System.err.println("Unsopported file: " + fileName);
			return TrackFormat.BED;
		}
	}
	
	public static LinkedHashMap<String, IntervalFeatureSet> createIntervalFeatureSets(List<String> fileNames) throws IOException{
		LinkedHashMap<String, IntervalFeatureSet> ifsets= new LinkedHashMap<String, IntervalFeatureSet>();
		for(String x : fileNames){
			File f= new File(x);
			if(getFileTypeFromName(x).equals(TrackFormat.BED) || getFileTypeFromName(x).equals(TrackFormat.GFF)){
				if(!ifsets.containsKey(x)){ // If the input has duplicates, do not reload duplicates!
					IntervalFeatureSet ifs= new IntervalFeatureSet(f);
					ifsets.put(x, ifs);
				}
			}
		}
		return ifsets;
	}
	
    /** 
     * Transpose list of list as if they were a table. No empty cells should be present. 
     * See http://stackoverflow.com/questions/2941997/how-to-transpose-listlist
     * FROM
     * [[a, b, c, d], [a, b, c, d], [a, b, c, d]] 
     * TO
     * [[a, a, a, a],
     *  [b, b, b, b],
     *  [c, c, c, c]]
     * @param table Table like list of lists, no empty cells.
     * @return
     */
    public static <T> List<List<T>> transpose(List<List<T>> table) {
    	List<List<T>> ret = new ArrayList<List<T>>();
        final int N = table.get(0).size();
        for (int i = 0; i < N; i++) {
            List<T> col = new ArrayList<T>();
            for (List<T> row : table) {
                col.add(row.get(i));
            }
            ret.add(col);
        }
        return ret;
    }
	
    /**
     * Map list of values mapping top genomic positions to a smaller list of positions by averaging 
     * values mapped to the same reference position. These averages are the values that will be used on the 
     * y-axis.
     * @param values Values to collapse
     * @param valuePositions Genomic position of the values
     * @param referencePositions Arrival positions. I.e. where the valuePositions should be mapped to. 
     * Typically this is obtained from GenomicCoords.getMapping();  
     * @return
     */
    public static List<Double> collapseValues(List<Double> values, List<Integer> valuePositions, 
    		List<Double> referencePositions){

    	// First store here all the values mapping to each reference position, then take the average.
    	LinkedHashMap<Integer, List<Double>> zlist= new LinkedHashMap<Integer, List<Double>>();
    	for(int i= 0; i < referencePositions.size(); i++){
    		zlist.put(i, new ArrayList<Double>());
    	} 
    	
    	for(int i= 0; i < valuePositions.size(); i++){
    		if(values.get(i) == null){
    			continue;
    		}
    		int pos= valuePositions.get(i);
    		if(pos >= referencePositions.get(0) && pos <= referencePositions.get(referencePositions.size()-1)){
	    	// Do not consider data points outside screenMap.
    			int j= Utils.getIndexOfclosestValue(pos, referencePositions);
    			zlist.get(j).add((double)values.get(i));
    		}
    	}
    	
    	List<Double> compressed= new ArrayList<Double>();
    	for(int i= 0; i < referencePositions.size(); i++){
    		compressed.add(Utils.calculateAverage(zlist.get(i)));
    	}
    	return compressed;
    }

	/**
	 * Get sequence as byte[] for the given genomic coords.
	 * @param fasta
	 * @param gc
	 * @return
	 * @throws IOException
	 */
	public static byte[] prepareRefSeq(String fasta, GenomicCoords gc) throws IOException{

		byte[] faSeq= null;
		if(fasta != null){
			IndexedFastaSequenceFile faSeqFile = null;
			try {
				faSeqFile = new IndexedFastaSequenceFile(new File(fasta));
				try{
					faSeq= faSeqFile.getSubsequenceAt(gc.getChrom(), gc.getFrom(), gc.getTo()).getBases();
				} catch (NullPointerException e){
					System.err.println("Cannot fetch sequence " + gc.toString());
					e.printStackTrace();
				}
				faSeqFile.close();
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		}
		return faSeq;
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
			throw new RuntimeException("Invalid string to convert to int: " + x);
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
		int halfWindow= (int)Math.rint(windowSize / 2d);
		if(rawInput.trim().equals("ff")){				
			from += halfWindow; 
			to += halfWindow;
			if(!gc.getSamSeqDict().isEmpty()){
				int chromLen= gc.getSamSeqDict().getSequence(chrom).getSequenceLength();
				if(to > chromLen){
					to= chromLen;
					from= to - gc.getUserWindowSize() + 1;
				}
			}			
			return chrom + ":" + from + "-" + to;
		} else if(rawInput.trim().equals("bb")) {
			from -= halfWindow;
			to -= halfWindow; 
			if(from < 1){
				from= 1;
				to= from + gc.getUserWindowSize() - 1;
			}
			return chrom + ":" + from + "-" + to;
		} else if(rawInput.trim().equals("f")){
			int step= (int)Math.rint(windowSize / 10d);

			from += step; 
			to += step;
			if(!gc.getSamSeqDict().isEmpty()){
				int chromLen= gc.getSamSeqDict().getSequence(chrom).getSequenceLength();
				if(to > chromLen){
					to= chromLen;
					from= to - gc.getUserWindowSize() + 1;
				}
			}			
			return chrom + ":" + from + "-" + to;
			//step= (step == 0) ? 1 : step;
			//from += step; 
			//to += step;
			//return chrom + ":" + from + "-" + to;
		} else if(rawInput.trim().equals("b")){
			int step= (int)Math.rint(windowSize / 10d);
			from -= step;
			to -= step; 
			if(from < 1){
				from= 1;
				to= from + gc.getUserWindowSize() - 1;
			}
			return chrom + ":" + from + "-" + to;
			//int step= (int)Math.rint(windowSize / 10d);
			//step= (step == 0) ? 1 : step;
			//from -= step; 
			//to -= step;
			//return chrom + ":" + from + "-" + to;
		} else if(rawInput.trim().matches("\\d+.*")) { // You might want to be more specific
			return chrom + ":" + parseGoToRegion(rawInput);
		} else if(rawInput.trim().startsWith("+") 
				|| rawInput.trim().startsWith("-")){
			int offset= parseStringToIntWithUnits(rawInput.trim());
			from += offset;
			if(from <= 0){
				from= 1;
				to= gc.getGenomicWindowSize();
			} else {
				to += offset;
			}
			return chrom + ":" + from + "-" + to;
		}else if (rawInput.equals("q")) {
			System.exit(0);	
		} else {
			throw new RuntimeException("Invalid input for " + rawInput);
		}
		return region;
	}
	
	/**Parse the rawInput string in the form ':123-456' to return either
	 * the first int or both ints. 
	 * */
	private static String parseGoToRegion(String rawInput){
		String[] fromTo= rawInput.trim().replaceAll(",", "").
								         replaceAll(" ", "").split("-");
		if(fromTo.length == 1){
			Integer.parseInt(fromTo[0]); // Check you actually got an int.
			return fromTo[0];
		} else {
			Integer.parseInt(fromTo[0]); // Check you actually got an int.
			Integer.parseInt(fromTo[1]);
			return fromTo[0] + "-" + fromTo[1];
		}
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
	
	public static List<Double> seqFromToLenOut(double from, double to, int lengthOut){
		
		if(lengthOut < 0){
			String msg= "Sequence from " + from + " to " + to + ". Invalid lenght of sequence: Cannot be < 1. Got " + lengthOut;
			throw new RuntimeException(msg);
		} 
		if(lengthOut == 0){ // Consistent with R seq(): return 0-length vector
			return new ArrayList<Double>();
		}		
		
		List<Double> mapping= new ArrayList<Double>();
		double span= to - from + 1;
		double step= ((double)span - 1)/(lengthOut - 1);
		mapping.add((double)from);
		for(int i= 1; i < lengthOut; i++){
			mapping.add((double)mapping.get(i-1)+step);
		}

		if(lengthOut == 1){ // Consistent with R seq(from, to, length.out= 1) -> from
		//	mapping.add((double)from);
			return mapping;
		}
		
		double diffTo= Math.abs(mapping.get(mapping.size() - 1) - to);
		if(diffTo > Math.abs((to + 1e-9))){
			String msg= "Error generating sequence from " + from + " to " + to + " length " + lengthOut + "\n" +
					     "Last point: " + mapping.get(mapping.size() - 1) + "\n" +
					     "To diff: " + diffTo + "\n" +
					     "Step: " + step;
			throw new RuntimeException(msg);
		} else {
			mapping.set(mapping.size()-1, (double)to);
		}
		
		double diffFrom= Math.abs(mapping.get(0) - from);		
		if(diffFrom > 0.01 || mapping.size() != lengthOut){
			String msg= "Error generating sequence:\n" +
					    "Expected size: " + lengthOut + "; Effective: " + mapping.size() + "\n" + 
					    "From diff: " + diffFrom;
			throw new RuntimeException(msg);
		}
		return mapping;
	}

	/*Nicely tabulate list of rows. Each row is tab separated **/
	public static List<String> tabulateList(List<String> rawList) {
		
		// * Split each row in a list of strings. I.e. make list of lists
		List<ArrayList<String>> rawTable= new ArrayList<ArrayList<String>>();
		int ncol= 0;
		for(String x : rawList){
			List<String> row = new ArrayList<String>();
			for(String item : x.split("\t")){
				row.add(item);
			}
			rawTable.add((ArrayList<String>) row);
			// * Get max number of columns (max list length)
			if(row.size() > ncol){
				ncol= row.size(); 
			}
		}
				
		// * Iterate through each column 
		List<ArrayList<String>> paddedTable= new ArrayList<ArrayList<String>>();
		
		for(int i= 0; i < ncol; i++){
			// Collect all items in column i in a list. I.e. select column i
			List<String> col= new ArrayList<String>();
			for(ArrayList<String> row : rawTable){
				if(row.size() > i){
					col.add(row.get(i));
				} else { // If a row doesn't have enough columns add a dummy field 
					col.add("");
				}
			}
			// Get the longest string in this column
			int maxStr= 0;
			for(String x : col){
				if(x.length() > maxStr){
					maxStr= x.length();
				}
			}
			// ** Pass thorugh the column again and pad with spaces to match length of longest string
			// maxStr+=1; // +1 is for separating
			for(int j= 0; j < col.size(); j++){
				String padded= String.format("%-" + maxStr + "s", col.get(j));
				col.set(j, padded);
			}
			paddedTable.add((ArrayList<String>) col);
		}
		// Each list in padded table is a column. We need to create rows as strings
		List<String> outputTable= new ArrayList<String>();
		for(int r= 0; r < paddedTable.get(0).size(); r++){
			StringBuilder row= new StringBuilder();
			for(int c= 0; c < paddedTable.size(); c++){
				row.append(paddedTable.get(c).get(r) + " ");
			}
			outputTable.add(row.toString().trim());
		}
		return outputTable;
	}		
}
