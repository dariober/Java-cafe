package tracks;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.lang3.StringUtils;

import com.google.common.base.Joiner;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import samTextViewer.Utils;
import samTextViewer.SamLocusIterator.LocusInfo;

/**
 * Create printable methylation track from list of LocusInfo and matched reference sequence.
 * @author berald01
 */
public class TrackMethylation {

	List<ScreenLocusInfo> screenLocusInfoList; 
	private double scorePerDot;
	
	/* C o n s t r u c t o r s */
	public TrackMethylation(List<ScreenLocusInfo> screenLocusInfoList){		
		this.screenLocusInfoList= screenLocusInfoList;
	}
	
	/* M e t h o d s */
	
	public String printToScreen(int yMaxLines, boolean noFormat){
		
		// Get list of mean M and U
		List<Double> mValues= new ArrayList<Double>();
		List<Double> uValues= new ArrayList<Double>();
		for(ScreenLocusInfo x : this.screenLocusInfoList){
			mValues.add(x.getMeanCntM());
			uValues.add(x.getMeanCntU());
		}		
		ArrayList<List<String>> profile= (ArrayList<List<String>>) getMethylProfileList(
				mValues, uValues, yMaxLines, noFormat);
		
		ArrayList<String> lineStrings= new ArrayList<String>();
		for(int i= (profile.size() - 1); i >= 0; i--){
			List<String> xl= profile.get(i);
			Set<String> unique= new HashSet<String>(xl);
			if(unique.size() == 1 && unique.contains(" ")){ // Do not print blank lines
				continue;
			} else {
				lineStrings.add(StringUtils.join(xl, ""));
			}
		}
		return Joiner.on("\n").join(lineStrings);
		
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
	* @param mValues, uValues Counts of methylated and unmethylated calls. MEMO: Type is double
	* since they might be averages of several positions mapping to the same screen position.
	 * @return
	 */
	private List<List<String>> getMethylProfileList(List<Double> mValues, List<Double> uValues, 
			int yMaxLines, boolean noFormat){
		
		if(mValues.size() != uValues.size()){ // Sanity check
			System.err.println("Number of M-values does not equal number of U-values");
			System.err.println("M-values: " + mValues);
			System.err.println("U-values: " + uValues);
			System.exit(1);
		}
		
		String charForM= (noFormat) ? "*" : "\033[31m*\033[0m"; // 107: white bg; 31: red text
		String charForU= (noFormat) ? "." : "\033[34m.\033[0m"; // 107: white bg; 34: blue text
		String charForZero= (noFormat) ? "_" : "\033[_\033[0m"; // 107: white bg
		String charForNonCyt= " ";
		String charForFill= " ";
		
		// Rescale counts to fit y-axis
		List<Double> rescaledM= new ArrayList<Double>();
		List<Double> rescaledU= new ArrayList<Double>();
		for(int i= 0; i < mValues.size(); i++){
			if(mValues.get(i) == null){
				rescaledM.add(0.0);
				rescaledU.add(0.0);
			} else {
				rescaledM.add(mValues.get(i) / this.getMaxDepth() * yMaxLines);
				rescaledU.add(uValues.get(i) / this.getMaxDepth() * yMaxLines);
			}
		}
		scorePerDot= this.getMaxDepth() / yMaxLines;
		List<List<String>> profile= new ArrayList<List<String>>();
		
		// Produce profile strings: Each loop is a vertical bar.
		for(int i= 0; i < rescaledM.size(); i++){
			ArrayList<String> strDepth= new ArrayList<String>(); // This will be a vertical bar
			double locDepth= rescaledM.get(i) + rescaledU.get(i);
			double locDepthM= rescaledM.get(i);
			double locDepthU= rescaledU.get(i);
			
			if(Double.isNaN(locDepth)){ // This position(s) are not a cytosines
				strDepth.add(charForNonCyt);
			} else if((int)Math.rint(locDepth) == 0){ // For zero coverage.
				strDepth.add(charForZero);
			} else {
				int nM= (int) Math.rint((locDepthM)); // How many M
				List<String> listM = Collections.nCopies(nM, charForM);
				int nU= (int) Math.rint((locDepthU)); // How many U
				List<String> listU = Collections.nCopies(nU, charForU);
				strDepth.addAll(listM);
				strDepth.addAll(listU);
			}
			while(strDepth.size() < yMaxLines){ // Fill up list with blanks
				strDepth.add(charForFill);
		    }
			if(strDepth.size() != yMaxLines){ // Sanity check: each vertical bar is yMaxLines high.
				System.err.println("Error producing methylation profile: " + strDepth);
				System.err.println("Input mValues: " + mValues);
				System.err.println("Input rescaled mValues: " + rescaledM);
				System.err.println("Input rescaled uValues: " + rescaledU);
				System.err.println("Input uValues: " + uValues);
				System.err.println("Expected height: " + yMaxLines + "; Got: " + locDepth);
				System.exit(1);
			}
			profile.add(strDepth);
		}
		return Utils.transpose(profile);
	}
	
	/**
	 * Max depth in the region, depth sum of M+U.
	 * @return
	 */
	public double getMaxDepth(){
		double maxDepth= 0;
		for(ScreenLocusInfo x : this.screenLocusInfoList){
			double depth= x.getMeanCntM() + x.getMeanCntU();
			if(maxDepth < depth){
				maxDepth= depth;
			}
		}
		return maxDepth;
	}

	public double getScorePerDot(){
		return scorePerDot;
	}
}

/*
 	public String printToScreen(List<Double> screenToGenomeMap, int yMaxLines, boolean noFormat){
		
		// Need to convert list of integer to double
		List<Double> mValues= new ArrayList<Double>();
		List<Double> uValues= new ArrayList<Double>();
		for(int i= 0; i < this.cntM.size(); i++){
			if(this.cntM.get(i) == null){
				mVal.add(null);
			} else {
				mVal.add((double)this.cntM.get(i));
			}
			if(this.cntU.get(i) == null){
				uVal.add(null);
			} else {
				uVal.add((double)this.cntU.get(i));
			}
		}			
		List<Double> mValues= Utils.collapseValues(mVal, this.positions, screenToGenomeMap);
		List<Double> uValues= Utils.collapseValues(uVal, this.positions, screenToGenomeMap);
		
		ArrayList<List<String>> profile= (ArrayList<List<String>>) getMethylProfileList(
				mValues, uValues, yMaxLines, noFormat);
		
		ArrayList<String> lineStrings= new ArrayList<String>();
		for(int i= (profile.size() - 1); i >= 0; i--){
			List<String> xl= profile.get(i);
			Set<String> unique= new HashSet<String>(xl);
			if(unique.size() == 1 && unique.contains(" ")){ // Do not print blank lines
				continue;
			} else {
				lineStrings.add(StringUtils.join(xl, ""));
			}
		}
		return Joiner.on("\n").join(lineStrings);
		
	}	

 	public TrackMethylation(List<LocusInfo> locusInfoList, IndexedFastaSequenceFile faSeqFile){

		for(LocusInfo loc : locusInfoList){
			positions.add(loc.getPosition());
			MethylLocus ml= new MethylLocus(loc, faSeqFile);
			cntM.add(ml.getCntM());
			cntU.add(ml.getCntU());
		}
	}

 * */