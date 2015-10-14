package methylationViewer;

// import htsjdk.samtools.util.SamLocusIterator.LocusInfo;
import samTextViewer.SamLocusIterator;
import samTextViewer.SamLocusIterator.LocusInfo;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Set;

import org.apache.commons.lang3.StringUtils;

import samTextViewer.Utils;
import coverageViewer.CoverageViewer;

/**
 * Class for working with methylated loci.
 * @author berald01
 *
 */
public class MethylLoci{
	
	private List<MethylLocus> methylLociList= new ArrayList<MethylLocus>();
	private String chrom; 
	private List<Integer> genomicPositions= new ArrayList<Integer>();
	private List<Float> MDepth= new ArrayList<Float>();
	private List<Float> UDepth= new ArrayList<Float>();
	private List<Character> refBases= new ArrayList<Character>();
	// private float bpPerChar; // Always upper case
	
	/*    C o n s t r u c t o r s    */
	public MethylLoci(List<MethylLocus> metlList){
		methylLociList= metlList;
		chrom= metlList.get(0).getLocus().getSequenceName();

		for(MethylLocus m : metlList){ // Check loci come from the same reference sequence
			if(m.getLocus().getSequenceName() != chrom){
				System.err.println("Mixing loci from different chroms is not allowed");
				System.exit(1);
			}
		}
		this.pupulateList();
	}

	/**
	 * Create list of methylLoci by extracting LocusInfo from CoverageViewer obj. 
	 * LocusInfo converted to MethylLocus using the given reference sequence and
	 * offset.  
	 * @param cw
	 * @param subRefSeq 
	 * @param offset subSeqRef starts at this genomic position. offset=1 if subRefSeq starts at the beginning
	 * of the chromosome.  
	 */
	public MethylLoci(CoverageViewer cw, byte[] subRefSeq, int offset){

		if(cw.getDepth().size() > subRefSeq.length){
			System.err.println("Reference sequence does not span the size of the view. View length: " + cw.getDepth().size() 
					+ " Ref length: " + subRefSeq.length);
			System.exit(1);
		}
		
		for(LocusInfo locus : cw.getLocusInfoList()){
			MethylLocus methylLocus= new MethylLocus(locus, subRefSeq, offset);
			methylLociList.add(methylLocus);
		}
		chrom= methylLociList.get(0).getLocus().getSequenceName();
		this.pupulateList();
	}
	
	
             	/*    M e t h o d s   */
	
	private void pupulateList(){
		
		for(MethylLocus ml : this.methylLociList){  

			this.genomicPositions.add(ml.getLocus().getPosition());
			this.refBases.add(Character.toUpperCase(ml.getRefBase()));
			if(ml.getCntM() == null && ml.getCntU() == null){
				this.MDepth.add(Float.NaN);
				this.UDepth.add(Float.NaN);
			} else {
				this.MDepth.add((float)ml.getCntM());
				this.UDepth.add((float)ml.getCntU());
			}

			if(genomicPositions.size() > 1){ // QC for loci in ascending order 
				if(genomicPositions.get(genomicPositions.size() - 1) <= genomicPositions.get(genomicPositions.size() - 2)){
					System.err.println("Methylated loci not in positional order!");
					System.exit(1);
				}
			}
		}
	}	

	/**
	 * List of strings where each string is a horizontal line of methylation track.
	 * Adapted from CoverageViewer.getProfileStrings(int ymaxLines)
	 * @param ymaxLines Rescale the y-axis (M+U) to be at most this many lines.
	 * @param noFormat apply format to string
	 * @return
	 */
	public List<String> getMethylProfileStrings(int ymaxLines, boolean noFormat){
		
		ArrayList<String> depthStrings= new ArrayList<String>();
		ArrayList<List<String>> profile= (ArrayList<List<String>>) 
				getMethylProfileList(this.MDepth, this.UDepth, ymaxLines, noFormat);
			
		for(int i= (profile.size() - 1); i >= 0; i--){
			List<String> xl= profile.get(i);
			Set<String> unique= new HashSet<String>(xl);
			if(unique.size() == 1 && unique.contains(CoverageViewer.FILL)){ // Do not print blank lines
				continue;
			} else {
				depthStrings.add(StringUtils.join(xl, ""));
			}
		}
		return depthStrings;
	}
	
	/**
	 * Produce a representation of methylation using text characters. 
	 * Output is a list of lists where each inner list is a horizontal line of 
	 * the methylation track. The vertical y-axis is scaled to be at most ymaxLines.
	 * Adapted from CoverageViewer.getProfileList(). It will look like
	 * [ 
	 *   ["M", "u", " ", "u", " "],
	 *   ["M", "u", " ", "u", " "],
	 *   ["M", "m", " ", "m", "U"],
	 *   ["M", "m", " ", "m", "U"],
	 *   ["M", "m", "m", "m", "M"],
	 * ]
	 * @param depth methylation along positions
	 * @param ymaxLines Max number of text lines to use. 
	 * Values in depthArray will be rescaled accordingly. 0 or -ve to disable scaling.
	 * @param noFormat Should the formatting be excluded?
	 * @return
	 */
	private List<List<String>> getMethylProfileList(List<Float> depthM, List<Float> depthU, 
			int ymaxLines, boolean noFormat){
		
		String charForM= (noFormat) ? "*" : "\033[107;31m*\033[0m"; // 107: white bg; 31: red text
		String charForU= (noFormat) ? "." : "\033[107;34m.\033[0m"; // 107: white bg; 34: blue text
		String charForZero= (noFormat) ? "_" : "\033[107m_\033[0m"; // 107: white bg
		if(ymaxLines > 0 && ymaxLines <= this.getMaxDepth()){ // Rescale depth as required
			for(int i= 0; i < depthM.size(); i++){
				float rescaledM= (float) depthM.get(i) / this.getMaxDepth() * ymaxLines;
				depthM.set(i, rescaledM); 
				float rescaledU= (float) depthU.get(i) / this.getMaxDepth() * ymaxLines;
				depthU.set(i, rescaledU); 
			}
		}
		List<List<String>> profile= new ArrayList<List<String>>();

		for(int i= 0; i < depthM.size(); i++){
			ArrayList<String> strDepth= new ArrayList<String>(); // This will be a vertical bar
			float locDepth= depthM.get(i) + depthU.get(i);
			float locDepthM= depthM.get(i);
			float locDepthU= depthU.get(i);
			
			if(Float.isNaN(locDepth)){
				strDepth.add(" ");
			} else if((int)Math.round(locDepth) == 0){ // For zero coverage.
				strDepth.add(charForZero);
			} else {
				int nM= Math.round((locDepthM)); // How many M
				List<String> listM = Collections.nCopies(nM, charForM);
				int nU= Math.round((locDepthU)); // How many U
				List<String> listU = Collections.nCopies(nU, charForU);
				strDepth.addAll(listM);
				strDepth.addAll(listU);
			}
			// Fill up list with blanks
			while(strDepth.size() < this.getMaxDepth()){
				strDepth.add(CoverageViewer.FILL);
		    }
			profile.add(strDepth);
		}
		return CoverageViewer.transpose(profile);
	}
	
	/**
	 * Max depth in the region, depth sum of M+U.
	 * @return
	 */
	public float getMaxDepth(){
		float maxDepth= 0;
		for(int i= 0; i < this.MDepth.size(); i++){
			float depth= this.MDepth.get(i) + this.UDepth.get(i);
			if(maxDepth < depth){
				maxDepth= depth;
			}
		}
		return maxDepth;
	}
	
    /**
     * Adapted from CoverageViewer.compressCovergeViewer()
     * Compress large interval to reduce it to number of elements given in "windowSize".
     * The new genomicPositions will correspond to the first position of each group of elements.  
     * @param windowSize Number of elements to reduce the track to. I.e. the number of chars to
     * display horizontally.
     */
    public void compressCovergeViewer(int windowSize){
		/* Non C/G sites have NaN counts of course and should be ignored for compression */
    	// * Compress MDepth
    	LinkedHashMap<Integer, Float> zMDepth= Utils.compressNumericListTO_BE_DEPRECATED(this.getMDepth(), windowSize);
    	// * Compress UDepth
    	LinkedHashMap<Integer, Float> zUDepth= Utils.compressNumericListTO_BE_DEPRECATED(this.getUDepth(), windowSize);
    	
    	this.MDepth= new ArrayList<Float>(zMDepth.values()); // Summarized values for each group. This will be the y-axis
		this.UDepth= new ArrayList<Float>(zUDepth.values());
		ArrayList<Integer> zidx= new ArrayList<Integer>(zMDepth.keySet()); // Indexes of the first element of each group.
		List<Integer> newDepthAt= new ArrayList<Integer>();
		//this..maxDepth= 0; // Recalculate max depth
		for(int i=0; i < zidx.size(); i++){
			newDepthAt.add(this.getGenomicPositions().get(zidx.get(i)));
		//	if(this.depth.get(i) > this.maxDepth){
		//		this.maxDepth= this.depth.get(i); 
		//	}
		}
		this.genomicPositions= newDepthAt;
		int to= this.getMethylLociList().get(this.getMethylLociList().size()-1).getLocus().getPosition();
		int from= this.getMethylLociList().get(0).getLocus().getPosition();
		// this.bpPerChar= (float)(to - from + 1) / zMDepth.size(); 
		
    }
	
	public String toString(){
		StringBuilder sb= new StringBuilder();
		sb.append(this.chrom + ":" + this.genomicPositions.get(0) + "-" +
				this.genomicPositions.get(this.genomicPositions.size()-1) + "\n");
		sb.append("N. loci: " + this.genomicPositions.size() + "\n");
		sb.append("Genomci positions: " + this.genomicPositions + "\n");
		sb.append("Reference: " + this.refBases + "\n");
		sb.append("Depth M: " + this.MDepth + "\n");
		sb.append("Depth U: " + this.UDepth + "\n");
		return sb.toString();
	}
	
	/*       S e t t e r s   and   G e t t e r s     */
	
	public List<MethylLocus> getMethylLociList() {
		return methylLociList;
	}

	public String getChrom() {
		return chrom;
	}

	public List<Integer> getGenomicPositions() {
		return genomicPositions;
	}

	public List<Float> getMDepth() {
		return MDepth;
	}

	public List<Float> getUDepth() {
		return UDepth;
	}

	public List<Character> getRefBases() {
		return refBases;
	}	
}