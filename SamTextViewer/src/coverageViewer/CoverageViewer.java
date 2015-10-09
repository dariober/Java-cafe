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
import samTextViewer.Utils;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.SamLocusIterator;
import htsjdk.samtools.util.SamLocusIterator.LocusInfo;

public class CoverageViewer {

	/** Max number of loci per window (column of text char). If you query 
	 * large intervals, only some loci will be stored. This const define the density of
	 * loci per window. */
	final static int LOC_PER_WINDOW= 100;
	//final static String DOT= "*";
	final static String FILL= " ";
	
	/** List of positions and corresponding depth. Each element need not to represent a single
	 * genomic bp. An element might be an summary (e.g. mean) of a group of adjacent positions. */
	// private List<Integer> depth= new ArrayList<Integer>();
	private List<Float> depth= new ArrayList<Float>();

	/** Max value in depth List */
	private float maxDepth= 0;
	
	/** Ruler: For each element in depth list say what position in genomic coordinates it corresponds to.
	 * Need not to be an ungapped sequence since each element in depth might be a group positions. In
	 * such case "depthAt" refers to the first base position of the group.*/
	private List<Integer> depthAt= new ArrayList<Integer>();
	
	private String chrom;
	private int from;
	private int to;
	/**
	 * Number of bases (bp) per character text on the screen. Set to 1 to start with
	 * Increased as large genomci windows are compressed.
	 * */
	private float bpPerChar= 1;
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
	public CoverageViewer(String sam, String chrom, int from, int to, int windowSize, 
			List<SamRecordFilter> filters){

		this.chrom= chrom;
		this.from= from;
		this.to= to;
		
		int range= to -from + 1;
		float density= 1/((float)range / ((float)windowSize * LOC_PER_WINDOW));
		
		SamReader samReader= ReadWriteBAMUtils.reader(sam, ValidationStringency.SILENT);
		SAMFileHeader fh= samReader.getFileHeader();
		
		IntervalList il= new IntervalList(fh);
		Random rand= new Random(); // Really you shouldn;t use rnd o filter loci.
		for(int i= from; i <= to; i++){
			float p= rand.nextFloat();
			if(p < density){
				Interval interval= new Interval(chrom, i, i);
				il.add(interval);
			}
		}
		SamLocusIterator samLocIter= new SamLocusIterator(samReader, il, true);
		samLocIter.setSamFilters(filters);
		Iterator<LocusInfo> iter= samLocIter.iterator();
	
		while(iter.hasNext()){
			LocusInfo locusInfo= iter.next();
			float curDepth= locusInfo.getRecordAndPositions().size();
			this.depth.add(curDepth);
			this.depthAt.add(locusInfo.getPosition());
			if(curDepth > this.maxDepth){
				this.maxDepth= curDepth;
			}
		}	
	}

	
	/**
	 * Initialize coverage track directly with lists of ints.
	 * @param depth
	 */
	public CoverageViewer(List<Float> depth, List<Integer> depthAt){
		this.depth= depth;
		this.depthAt= depthAt;
		float maxDepth= 0;
		for(float x : depth){
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
	 * @param depth Depth along positions
	 * @param ymaxLines Max number of text lines to use. 
	 * Values in depthArray will be rescaled accordingly. 0 or -ve to disable scaling.
	 * @return
	 */
	private List<List<String>> getProfileList(List<Float> depth, int ymaxLines){
		
		if(ymaxLines > 0 && ymaxLines <= this.maxDepth){ // Rescale depth as required
			for(int i= 0; i < depth.size(); i++){
				float rescaled= (float) depth.get(i) / maxDepth * ymaxLines;
				depth.set(i, rescaled);  
			}
		}
		//float scailingFactor= (float)ymaxLines/((float)maxDepth + (float)0.0001);
		//for(int i= 0; i < depth.size(); i++){
		//	float rescaled= depth.get(i) * scailingFactor;
		//	depth.set(i, rescaled);  
		//}
		List<List<String>> profile= new ArrayList<List<String>>();
		
		for(int i= 0; i < depth.size(); i++){
			ArrayList<String> strDepth= new ArrayList<String>(); // This will be a vertical bar
			float locDepth= depth.get(i);
			
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
    private static <T> List<List<T>> transpose(List<List<T>> table) {
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
     * Compress coverage viewer track to reduce it to number of elements given in "windowSize".
     * The new depthAt positions will correspond to the first position of each group of elements.  
     * @param windowSize Number of elements to reduce the track to. I.e. the number of chars to
     * display horizontally.
     */
    public void compressCovergeViewer(int windowSize){
		LinkedHashMap<Integer, Float> zcw= Utils.compressListOfInts(this.getDepth(), windowSize);
		this.depth= new ArrayList<Float>(zcw.values()); // Summarized values for each group. This will be the y-axis
		ArrayList<Integer> zidx= new ArrayList<Integer>(zcw.keySet()); // Indexes of the first element of each group.
		List<Integer> newDepthAt= new ArrayList<Integer>();
		this.maxDepth= 0; // Recalculate max depth
		for(int i=0; i < zidx.size(); i++){
			newDepthAt.add(this.getDepthAt().get(zidx.get(i)));
			if(this.depth.get(i) > this.maxDepth){
				this.maxDepth= this.depth.get(i); 
			}
		}
		this.bpPerChar= (float)(this.to - this.from + 1) / zcw.size(); 
		this.depthAt= newDepthAt;
    }
        
    public String ruler(int markDist){
    	String numberLine= "";
    	int prevLen= 0;
    	int i= 0;
		while(i < this.depthAt.size()){
			String posMark= String.valueOf(this.getDepthAt().get(i));
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
    }
    
    public String toString(){
    	StringBuilder sb= new StringBuilder();
    	sb.append("depth: " + this.depth + "\n");
    	sb.append("depthAt: " + this.depthAt + "\n");
    	sb.append("maxDepth: " + this.maxDepth + "\n");
    	return sb.toString();
    }

    /* Setters and setters */
    
	public List<Float> getDepth(){
		return this.depth;
	}
    
	public float getMaxDepth() {
		return maxDepth;
	}

	public List<Integer> getDepthAt() {
		return depthAt;
	}

	public void setDepthAt(List<Integer> depthAt) {
		this.depthAt = depthAt;
	}

	public float getBpPerChar() {
		return bpPerChar;
	}

	public void setBpPerChar(float bpPerChar) {
		this.bpPerChar = bpPerChar;
	}
}
