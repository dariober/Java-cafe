package samTextViewer;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

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
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;
import java.util.zip.GZIPInputStream;

import org.apache.commons.lang3.text.StrTokenizer;
import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.tdf.TDFReader;

import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import tracks.IntervalFeatureSet;
import tracks.Track;
import tracks.TrackReads;

/**
 * @author berald01
 *
 */
public class Utils {
	
	/** Get the first chrom string from first line of input file. As you add support for more filetypes you should update 
	 * this function. This method is very dirty and shouldn't be trusted 100% */
	public static String initRegionFromFile(String x) throws IOException{
		String region= "";
		if(x.toLowerCase().endsWith(".bam") || x.toLowerCase().endsWith(".cram")){
			SamReaderFactory srf=SamReaderFactory.make();
			srf.validationStringency(ValidationStringency.SILENT);
			SamReader samReader = srf.open(new File(x));
			region= samReader.getFileHeader().getSequence(0).getSequenceName();
			samReader.close();
			return region;
		} else if(x.toLowerCase().endsWith(".bigwig") || x.toLowerCase().endsWith(".bw")){
			BBFileReader reader= new BBFileReader(x);
			region= reader.getChromosomeNames().get(0);
			reader.close();
			return region;
		} else if(x.toLowerCase().endsWith(".tdf")){
			Iterator<String> iter = TDFReader.getReader(x).getChromosomeNames().iterator();
			while(iter.hasNext()){
				region= iter.next();
				if(!region.equals("All")){
					return region;
				}
			} 
			System.err.println("Cannot initialize from " + x);
			throw new RuntimeException();
		} else if(x.toLowerCase().endsWith(".gz")){
			InputStream fileStream = new FileInputStream(x);
			GZIPInputStream gzipStream = new GZIPInputStream(fileStream);
			Reader decoder = new InputStreamReader(gzipStream, "UTF-8");
			BufferedReader br = new BufferedReader(decoder);
			while(true){
				region= br.readLine().trim();
				if(region == null){
					break;
				}
				region= (region.split("\t")[0]).trim();
				if(region.startsWith("#") || region.isEmpty()){
					continue;
				} else {
					br.close();
					return region;
				}
			}
			br.close();
			return region;
		} else {
			BufferedReader br = new BufferedReader(new FileReader(x));
			while(true){
				region= br.readLine().trim();
				if(region == null){
					break;
				}
				region= (region.split("\t")[0]).trim();
				if(region.startsWith("#") || region.isEmpty()){
					continue;
				} else {
					br.close();
					return region;
				}
			}
			br.close();
			return region;
		}
	}
	
	public static boolean bamHasIndex(String bam) throws IOException{

		SamReaderFactory srf=SamReaderFactory.make();
		srf.validationStringency(ValidationStringency.SILENT);
		SamReader samReader = srf.open(new File(bam));
		boolean hasIndex= samReader.hasIndex();
		samReader.close();
		return hasIndex;
		
	}
	
	public static GenomicCoords findNextStringOnFile(String string, String filename, GenomicCoords curGc,
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
	}
	
	public static GenomicCoords goToNextFeatureOnFile(String filename, GenomicCoords curGc, 
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
		return intervalFiles.get(chosenFn).coordsOfNextFeature(curGc);
	}
	 
	
	public static String getFileTypeFromName(String fileName){
		fileName= fileName.toLowerCase();
		
		if(    fileName.endsWith(".bed") 
		    || fileName.endsWith(".bed.gz") 
		    || fileName.endsWith(".bed.gz.tbi")){
			return "bed";
		} else if( fileName.endsWith(".gtf") 
				|| fileName.endsWith(".gtf.gz")
				|| fileName.endsWith(".gtf.gz.tbi")
				|| fileName.endsWith(".gff") 
				|| fileName.endsWith(".gff.gz") 
				|| fileName.endsWith(".gff.gz.tbi")){
			return "gff";
		} else if(fileName.endsWith(".bam") || fileName.endsWith(".cram")){
			return "bam";
		} else if(fileName.endsWith(".bigwig") || fileName.endsWith(".bw")) {
			return "bigWig";
		} else if(fileName.endsWith(".tdf")) {
			return "tdf";
		} else if(fileName.endsWith(".bedgraph.gz") || fileName.endsWith(".bedgraph")) {
			return "bedGraph";
		} else {
			// System.err.println("Unsopported file: " + fileName);
			return "bed";
		}
	}
	
	public static LinkedHashMap<String, IntervalFeatureSet> createIntervalFeatureSets(List<String> fileNames) throws IOException{
		LinkedHashMap<String, IntervalFeatureSet> ifsets= new LinkedHashMap<String, IntervalFeatureSet>();
		for(String x : fileNames){
			File f= new File(x);
			if(getFileTypeFromName(x).equals("bed") || getFileTypeFromName(x).equals("gff")){
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
	
	/** Get the coordinates of the first mapped read in bam file and put it
	 * in GenomicCoords obj.
	 * @param bam
	 */
	/*public static GenomicCoords getStartCoordsOfBAM(String bam){
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
	}*/
	/*public static GenomicCoords getStartCoordsOfBAM(String bam, String chrom){
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
	}*/	
	
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
		} else if(rawInput.trim().startsWith(":")) { // You might want to be more specific than just startsWith(:)
			String pos= rawInput.trim().replaceFirst(":", "");
			Integer.parseInt(pos); // Check you actually got an int.
			return chrom + ":" + pos;
		} else if(rawInput.trim().startsWith("+") 
				|| rawInput.trim().startsWith("-") 
				|| Character.isDigit(rawInput.trim().charAt(0))){
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
		
		if(lengthOut < 1){
			String msg= "Invalid lenght of sequence: Cannot be < 1. Got " + lengthOut;
			throw new RuntimeException(msg);
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
	
	/** From cmdInput extract regex and ylimits then iterate through the tracks list to set 
	 * the ylimits in the tracks whose filename matches the regex.
	 * The input list is updated in place! 
	*/
	public static List<Track> setTrackYlimitsForRegex(String cmdInput, List<Track> tracks) throws InvalidCommandLineException{

		StrTokenizer str= new StrTokenizer(cmdInput);
		str.setQuoteChar('\'');
		List<String> tokens= str.getTokenList();
		if(tokens.size() != 4){
			System.err.println("Error in :ylim subcommand. Expected 4 args got: " + cmdInput);
			throw new InvalidCommandLineException();
		}
		String ylimRegex= tokens.get(0);
		
		try{
			Pattern.compile(ylimRegex); // Validate regex
		} catch(PatternSyntaxException e){
	    	System.err.println("Invalid regex in: " + cmdInput);
	    	System.err.println(e.getDescription());
		}
		
		double ymin= Double.NaN;
        double ymax= Double.NaN;
		try{
			ymin= Double.parseDouble(tokens.get(1));
			ymax= Double.parseDouble(tokens.get(2));
		} catch(NumberFormatException e){
			ymin= Double.NaN;
			ymax= Double.NaN;
		}
		if(ymin >= ymax){
			System.err.println("Warning ymin >= ymax. Resetting to default.");
			ymin= Double.NaN;
			ymax= Double.NaN;							
		}
		for(Track tr : tracks){
			if(tr.getFilename().matches(ylimRegex)){
				tr.setYmin(ymin);
				tr.setYmax(ymax);
				//if(!(tr instanceof TrackReads)){
				//	String title= tr.getTitle().trim() + " ymin: " + tr.getYmin() + "; ymax: " + tr.getYmax() + ";\n";
				//	tr.setTitle(title);
				//}
			}
		}
		return tracks;
	}
	
	/**
	 * Generate sequence of doubles of desired length. Same as R seq(from, to, length.out)
	 * @param from
	 * @param to
	 * @param lengthOut
	 * @return
	 */
	/* public static List<Double> seqFromToLenOut(int from, int to, int lengthOut){
		
		if(lengthOut < 1){
			String msg= "Invalid lenght of sequence: Cannot be < 1. Got " + lengthOut;
			throw new RuntimeException(msg);
		}
		
		List<Double> mapping= new ArrayList<Double>();
		
		if(lengthOut == 1){ // Consistent with R seq(from, to, length.out= 1) -> from
			mapping.add((double)from);
			return mapping;
		}
		
		int span= to - from + 1;
		double step= ((double)span - 1)/(lengthOut - 1);
		mapping.add((double)from);
		for(int i= 1; i < lengthOut; i++){
			mapping.add((double)mapping.get(i-1)+step);
		}
		
		// First check last point is close enough to expection. If so, replace last point with
		// exact desired.
		double diffTo= Math.abs(mapping.get(mapping.size() - 1) - to);
		if(diffTo > (to * 0.001)){
			String msg= "Error generating sequence from " + from + " to " + to + " with length " + lengthOut + "\n" + 
					     "Last point: " + mapping.get(mapping.size() - 1) + "\n" +
					     "To diff: " + diffTo + "\n" +
					     "Step: " + step + "\n"
					     + "Sequence: " + mapping;
			
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
	} */
		
}
