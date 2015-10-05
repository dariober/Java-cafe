package coverageViewer;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.apache.commons.lang3.StringUtils;

import readWriteBAMUtils.ReadWriteBAMUtils;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.SamLocusIterator;
import htsjdk.samtools.util.SamLocusIterator.LocusInfo;

public class CoverageViewer {

	final static String DOT= "*";
	final static String FILL= " ";
	
	private List<Integer> depth;
//	private List<String> depthStrings;
	private int maxDepth= 0; // Max depth in interval 
//	private float depthPerLine; // After scaling, 1 line corresponds to 'depthPerLine' reads. 
	
	/* C o n s t r u c t o r s */

	public CoverageViewer(){
		
	}
	
	/**
	 * Construct coverage track from bam file and coordinates.
	 * @param sam
	 * @param chrom
	 * @param from
	 * @param to
	 */
	public CoverageViewer(String sam, String chrom, int from, int to){

		SamReader samReader= ReadWriteBAMUtils.reader(sam, ValidationStringency.SILENT);
		SAMFileHeader fh= samReader.getFileHeader();
		
		IntervalList il= new IntervalList(fh);
		Interval interval= new Interval(chrom, from, to);
		il.add(interval);
		SamLocusIterator samLocIter= new SamLocusIterator(samReader, il, true);
		
		// TODO: Implement SamRecordFilter by interpreting -f <int> and -F <int> flags
		
		List<Integer> depth= new ArrayList<Integer>();
		Iterator<LocusInfo> iter= samLocIter.iterator();		
		while(iter.hasNext()){
			int curDepth= iter.next().getRecordAndPositions().size();
			depth.add(curDepth);
			if(curDepth > this.maxDepth){
				this.maxDepth= curDepth;
			}
		}
		this.depth= depth;
		
	}
	
	/**
	 * Initialize coverage track directly with list of ints.
	 * @param depth
	 */
	public CoverageViewer(List<Integer> depth){
		this.depth= depth;
		int maxDepth= 0;
		for(int x : depth){
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
			if(unique.size() == 1 && unique.contains(FILL)){ // Do not print blank lines made of blanks.
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
	private List<List<String>> getProfileList(List<Integer> depth, int ymaxLines){
		
		if(ymaxLines > 0 && ymaxLines <= this.maxDepth){ // Rescale depth as required
			for(int i= 0; i < depth.size(); i++){
				int rescaled= (int) Math.round( (double) depth.get(i) / maxDepth * ymaxLines);
				depth.set(i, rescaled);  
			}
		}
		List<List<String>> profile= new ArrayList<List<String>>();
		
		for(int i= 0; i < depth.size(); i++){
			ArrayList<String> strDepth= new ArrayList<String>();
			int locDepth= depth.get(i);
			for(int j= 0; j < locDepth; j++){
				strDepth.add(DOT);
			}
			// Fill up list with blanks
			while(strDepth.size() < maxDepth){
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


    /**STUB: USed to filter out records. Doesn't seem to be very friendly though. 
     * Maybe better to write out a tmp bam with the required records?
     * @param f_incl
     * @param F_excl
     * @return
     */
    private List<SamRecordFilter> bitflagToSamRecFilterList(int f_incl, int F_excl){
		
    	final int[] bits= {2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048};
    	for(int i : bits){
    		if((f_incl & 2) == 2){
//    			new FailsVendorReadQualityFilter();
    		}
    	}
    	return null;
    }

    /* Setters and setters */
    
	public List<Integer> getDepth(){
		return this.depth;
	}
    
	public int getMaxDepth() {
		return maxDepth;
	}

}

/*
private void setProfileStrings(int ymaxLines){
	
	ArrayList<String> depthStrings= new ArrayList<String>();
	ArrayList<List<String>> profile= (ArrayList<List<String>>) getProfileList(this.depth, ymaxLines);
		
	for(int i= (profile.size() - 1); i >= 0; i--){
		List<String> xl= profile.get(i);
		Set<String> unique= new HashSet<String>(xl);
		if(unique.size() == 1 && unique.contains(FILL)){ // Do not print blank only lines
			continue;
		} else {
			depthStrings.add(StringUtils.join(xl, ""));
		}
	}
	this.depthStrings= depthStrings;
} */
