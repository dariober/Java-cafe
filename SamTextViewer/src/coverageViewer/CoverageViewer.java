package coverageViewer;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Random;
import java.util.Set;

import org.apache.commons.lang3.StringUtils;

import readWriteBAMUtils.ReadWriteBAMUtils;
import samTextViewer.GenomicCoords;
import samTextViewer.SamLocusIterator;
import samTextViewer.SamLocusIterator.LocusInfo;
import samTextViewer.Utils;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;

/**
 * This class should hold the information to display coverage and methylation profile
 * for given genomic coordinates.
 * 
 * @author berald01
 *
 */
public class CoverageViewer {

	/** Max number of loci per window (column of text char). If you query 
	 * large intervals, only some loci will be stored. This const define the density of
	 * loci per window. */
	final static int LOC_PER_WINDOW= 200;

	final static public String FILL= " ";
	
	/** Max value in depth List */
	private Double maxDepth= 0.0;
	
	/** For each element in depth list say what position in genomic coordinates it corresponds to. */
	//private List<Double> depth= new ArrayList<Double>();
	//private List<Integer> genomicPositions= new ArrayList<Integer>();

	private LinkedHashMap<Integer, Integer> depthAtPosition= new LinkedHashMap<Integer, Integer>();
	
	private List<Double> mappingToScreen;
	
	ArrayList<LocusInfo> locusInfoList= new ArrayList<LocusInfo>();; // This is going to be useful for methylation calling 
	
	/* C o n s t r u c t o r s */

	public CoverageViewer(){
		
	}
	
	/**
	 * Construct coverage track from bam file and coordinates.
	 * @param sam
	 * @param chrom
	 * @param from
	 * @param to
	 * @param windowSize
	 * @param filters
	 */
	public CoverageViewer(String sam, GenomicCoords gc, int windowSize,
			List<SamRecordFilter> filters){

		int range= gc.getTo() - gc.getFrom() + 1;
		double density= 1/((double)range / ((double)windowSize * LOC_PER_WINDOW));
		
		SamReader samReader= ReadWriteBAMUtils.reader(sam, ValidationStringency.SILENT);
		SAMFileHeader fh= samReader.getFileHeader();
		
		IntervalList il= new IntervalList(fh);
		Random rand= new Random(); // Really you shouldn;t use rnd o filter loci.
		for(int i= gc.getFrom(); i <= gc.getTo(); i++){
			float p= rand.nextFloat();
			if(p < density){
				Interval interval= new Interval(gc.getChrom(), i, i);
				il.add(interval);
			}
		}
		SamLocusIterator samLocIter= new SamLocusIterator(samReader, il, true);
		samLocIter.setSamFilters(filters);
		Iterator<samTextViewer.SamLocusIterator.LocusInfo> iter= samLocIter.iterator();
	
		while(iter.hasNext()){
			samTextViewer.SamLocusIterator.LocusInfo locusInfo= iter.next();
			locusInfoList.add(locusInfo);
			int curDepth= locusInfo.getRecordAndPositions().size();
			depthAtPosition.put(locusInfo.getPosition(), curDepth);
			//this.depth.add((double) curDepth);
			//this.genomicPositions.add(locusInfo.getPosition());
			if(curDepth > this.maxDepth){
				this.maxDepth= (double) curDepth;
			}
		}	
	}

	/**
	 * Initialize coverage track directly with list.
	 * @param depth
	 */
	public CoverageViewer(List<Double> depth, List<Integer> genomicPositions){
		this.depth= depth;
		this.genomicPositions= genomicPositions;
		Double maxDepth= 0.0;
		for(Double x : depth){
			if(x > maxDepth){
				maxDepth= x;
			}
		}
		this.maxDepth= maxDepth;
	}
	
	/* Methods */

	/**
	 * List of strings where each inner string is a horizontal line of coverage track.
	 * @param ymaxLines Rescale the y-axis (depth) to be at most this many lines.
	 * @return
	 */
	public List<String> getProfileStrings(int ymaxLines){
		
		// TODO? If the printable list is made only of zeros expand sites with decent number of
		// counts. E.g. outlier windows don't get squashed to zero. 
		// You could make the top 5% windows to be at least a . or a :
		
		ArrayList<String> depthStrings= new ArrayList<String>();
		ArrayList<List<String>> profile= (ArrayList<List<String>>) getProfileList(this.depth, ymaxLines);
			
		for(int i= (profile.size() - 1); i >= 0; i--){
			List<String> xl= profile.get(i);
			Set<String> unique= new HashSet<String>(xl);
			if(unique.size() == 1 && unique.contains(FILL)){ // Do not print blank lines
				continue;
			} else {
				depthStrings.add(StringUtils.join(xl, ""));
			}
		}
		return depthStrings;
	}
	

	/**
	 * Produce a representation of depth using text characters. 
	 * Output is a list of lists where each inner list is a horizontal line of 
	 * the coverage track. The vertical y-axis is scaled to be at most ymaxLines.
	 * @param depth2 Depth along positions
	 * @param ymaxLines Max number of text lines to use. 
	 * Values in depthArray will be rescaled accordingly. 0 or -ve to disable scaling.
	 * @return
	 */
	private List<List<String>> getProfileList(List<Double> depth2, int ymaxLines){
		
		ymaxLines= ymaxLines * 2; // Since we use : for 2x in a single line.
		
		if(ymaxLines > 0 && ymaxLines <= this.maxDepth){ // Rescale depth as required
			for(int i= 0; i < depth2.size(); i++){
				double rescaled= (double) depth2.get(i) / maxDepth * ymaxLines;
				depth2.set(i, rescaled);  
			}
		}
		List<List<String>> profile= new ArrayList<List<String>>();
		
		for(int i= 0; i < depth2.size(); i++){
			ArrayList<String> strDepth= new ArrayList<String>(); // This will be a vertical bar
			double locDepth= depth2.get(i);
			
			if((int)Math.round(locDepth) == 0){ // For zero coverage
				strDepth.add("_"); 
			}
			int nDouble= ((int)locDepth) / 2; // how many :
			for(int j= 0; j < nDouble; j++){
				strDepth.add(":");
			} 
			int diff= (int)Math.round(locDepth- (nDouble*2)); // Get remainder
			if(diff == 2){
				strDepth.add(":");
			} else if(diff == 1) {
				strDepth.add(".");
			} else if(diff == 0){
				//
			} else {
				System.err.println("Unexpected division");
				System.exit(1);
			}
			// Fill up list with blanks
			while(strDepth.size() < this.maxDepth){
				strDepth.add(FILL);
		    }
			profile.add(strDepth);
		}
		return transpose(profile);
	}
	
    /** 
     * Transpose list of list as if they were a table. No empty cells should be 
     * present. 
     * See http://stackoverflow.com/questions/2941997/how-to-transpose-listlist
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
     * Map flats in longList to the closest in referenceList. Elements mapping to the target
     * are summarised by averaging or other stats.
     * @param longList
     * @param referenceList
     * @return
     */
    public void setMappingToScreen(List<Double> screenMap){

    	LinkedHashMap<Integer, List<Double>> zlist= new LinkedHashMap<Integer, List<Double>>();
    	for(int i= 0; i < screenMap.size(); i++){
    		zlist.put(i, new ArrayList<Double>());
    	} 
    	
    	for(int i= 0; i < this.genomicPositions.size(); i++){
    		int pos= this.genomicPositions.get(i);
    		if(pos >= screenMap.get(0) && pos <= screenMap.get(screenMap.size()-1)){
	    	// Do not consider data points outside screenMap.
    			int j= Utils.getIndexOfclosestValue(pos, screenMap);
    			zlist.get(j).add(this.depth.get(i));
    		}
    	}
    	
    	List<Double> compressed= new ArrayList<Double>();
    	for(int i= 0; i < screenMap.size(); i++){
    		compressed.add(Utils.calculateAverage(zlist.get(i)));
    	}
    	this.mappingToScreen= compressed;
    }

    public List<Double> getMappingToScreen(List<Double> screenMap){
    	this.setMappingToScreen(screenMap);
    	return this.mappingToScreen;
    }
    
    /**
     * Compress coverage viewer track to reduce it to number of elements given in "windowSize".
     * The new genomicPositions positions will correspond to the first position of each group of elements.  
     * @param windowSize Number of elements to reduce the track to. I.e. the number of chars to
     * display horizontally.
     */
   /* public void compressCovergeViewer(int windowSize){
		LinkedHashMap<Integer, Float> zcw= Utils.compressNumericListTO_BE_DEPRECATED(this.getDepth(), windowSize);
		this.depth= new ArrayList<Float>(zcw.values()); // Summarized values for each group. This will be the y-axis
		ArrayList<Integer> zidx= new ArrayList<Integer>(zcw.keySet()); // Indexes of the first element of each group.
		List<Integer> newDepthAt= new ArrayList<Integer>();
		this.maxDepth= 0; // Recalculate max depth
		for(int i=0; i < zidx.size(); i++){
			newDepthAt.add(this.getGenomicPositions().get(zidx.get(i)));
			if(this.depth.get(i) > this.maxDepth){
				this.maxDepth= this.depth.get(i); 
			}
		}
		this.bpPerChar= (float)(this.to - this.from + 1) / zcw.size(); 
		this.genomicPositions= newDepthAt;
    } */
        
    /* public String ruler(int markDist){
    	String numberLine= "";
    	int prevLen= 0;
    	int i= 0;
		while(i < this.genomicPositions.size()){
			String posMark= String.valueOf(this.getGenomicPositions().get(i));
			if(i == 0){
				numberLine= posMark;
				i += posMark.length();
			} else if((numberLine.length() - prevLen) >= markDist){
				prevLen= numberLine.length();
				numberLine= numberLine + posMark;
				i += posMark.length();
			} else {
				numberLine= numberLine + " ";
				i++;
			}
		}
    	return numberLine;
    } */
    
    public String toString(){
    	StringBuilder sb= new StringBuilder();
    	sb.append("depth: " + this.depth + "\n");
    	sb.append("genomicPositions: " + this.genomicPositions + "\n");
    	sb.append("maxDepth: " + this.maxDepth + "\n");
    	return sb.toString();
    }

    /* Setters and setters */
    
	public List<Double> getDepth(){
		return this.depth;
	}
    
	public double getMaxDepth() {
		return maxDepth;
	}

	public List<Integer> getGenomicPositions() {
		return genomicPositions;
	}

	public void setGenomicPositions(List<Integer> depthAt) {
		this.genomicPositions = depthAt;
	}

	public ArrayList<LocusInfo> getLocusInfoList() {
		return locusInfoList;
	}
}
