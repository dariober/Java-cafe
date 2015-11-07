package tracks;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import samTextViewer.Utils;

/** Text representation of a continuous profile along the screen positions */
class TextProfile {
	
	private double maxDepth;   // Store the max depth of the track
	private double scorePerDot; // Store the scaling factor: Each dot in the profile cooresponds to
	                            // this many units of yValues. 
	// Text representation of yValues scled by yMaxLines. Each inner list is a line on screen.
	private List<List<String>> profile= new ArrayList<List<String>>(); 
	
	private String strFor1u= ".";
	private String strFor1uRev= "'";
	private String strFor1stNegU= ",";
	private String strFor2u= ":";
	private String strForFill= " ";
	private String strForZero= "_";
	private String strForZeroTop= " "; // Character.toString ((char) 773); // Upperscore, opposite of _
	// private String strForNaN= " ";
	
	/* C o n s t r u c t o r */
	/**
	 * @param yValues Values on the y-axis. Should be of the same length as the windowSize
	 * @param yMaxLines The yValues will be rescaled to fit this many lines of text.
	 */
	public TextProfile(List<Double> yValues, int yMaxLines){
		
		this.maxDepth= Integer.MIN_VALUE;;
        double minDepth= Integer.MAX_VALUE; 
		for(double y : yValues){
			if(y > this.maxDepth){
				this.maxDepth= y;
			}
			if(y < minDepth){
				minDepth= y;
			}
		}
		if(minDepth > 0){
			minDepth= 0;
		}
		if(maxDepth < 0){
		 	maxDepth= 0;
		}

		this.scorePerDot= (double)(this.maxDepth - minDepth) / (yMaxLines * 2); // * 2 because we use ':' for 2 units in a single line.
		
		// Locate zero on y axis. It's silly to generate a sequence just to find the index closest to zero. But anyway...
		List<Double> yAxis = Utils.seqFromToLenOut(minDepth, maxDepth, yMaxLines);
		int y0= Utils.getIndexOfclosestValue(0, yAxis);

		List<List<String>> profile= new ArrayList<List<String>>();
		for(int i= 0; i < yValues.size(); i++){
			double y= yValues.get(i);
			List<String> strDepth = prepareYColumn(y, yMaxLines, y0, this.scorePerDot);
			profile.add(strDepth);
			yColumnSanityCheck(strDepth, y0, y);
		}
		this.profile= Utils.transpose(profile);
	}

	/** Prepare a list of strings representing vertical bar. Bar height is yValue, rescaled to fit a y span of 
	 * yMaxLines of text. 
	 * @param yValue
	 * @param y0 index position of 0.
	 * @param scorePerDot
	 * @return
	 */
	private List<String> prepareYColumn(double yValue, int yMaxLines, int y0, double scorePerDot){

		Double yPosDotU= Math.abs(yValue / scorePerDot); // Y positions in dot units, not line units. 
		if((int)Math.rint(yPosDotU) == 0){ // For zero coverage
			ArrayList<String> strDepth= new ArrayList<String>();
			for(int j= 0; j < yMaxLines; j++){
				strDepth.add(strForFill);
			}			
			if(y0 < yMaxLines-1){
				strDepth.set(y0, this.strForZero);
			} else {
				strDepth.set(y0, this.strForZeroTop);
			}
			return strDepth;
		} else if(yValue < 0){
			return pileForNegative(yValue, yMaxLines, y0, scorePerDot);
		} else if(yValue > 0){
			return pileForPositive(yValue, yMaxLines, y0, scorePerDot);
		} else {
			throw new RuntimeException("Unexpected exception");
		}
	}
	
	private List<String> pileForPositive(double yValue, int yMaxLines, int y0, double scorePerDot){
		
		ArrayList<String> strDepth= new ArrayList<String>(); // This will be a vertical bar
		for(int j= 0; j < yMaxLines; j++){
			strDepth.add(this.strForFill);
		}
		int pos= y0;
		double currentScore= 0;
		while(true){
			if((yValue - currentScore) > scorePerDot * 1.5){ // Add double
				strDepth.set(pos, this.strFor2u);
				currentScore += 2 * scorePerDot;
			} else if((yValue - currentScore) > scorePerDot * 0.5){
				strDepth.set(pos, this.strFor1u);
				break;
			} else {
				break;
			}
			pos++;
		}
		return strDepth;
	}
	
	private List<String> pileForNegative(double yValue, int yMaxLines, int y0, double scorePerDot){
		
		ArrayList<String> strDepth= new ArrayList<String>(); // This will be a vertical bar
		for(int j= 0; j < yMaxLines; j++){
			strDepth.add(this.strForFill);
		}
		int pos= y0;
		double currentScore= 0;
		if(y0 < strDepth.size()-1){ // First char is 1u, unless all the values are negative and the zero is on the top
			strDepth.set(y0, this.strFor1stNegU);
			currentScore= 1 * scorePerDot;
			pos= y0-1;
		}
		while(pos >= 0){
			
			if((-yValue - currentScore) > scorePerDot * 1.5){ // Add double
				strDepth.set(pos, this.strFor2u);
				currentScore += 2 * scorePerDot;
			} else if((-yValue - currentScore) > scorePerDot * 0.5){
				strDepth.set(pos, this.strFor1uRev);
				break;
			} else {
				break;
			}
			pos--;
		}
		return strDepth;
	}
	
	private void yColumnSanityCheck(List<String> strDepth, int y0, double yValue){

		if(yValue < 0){ // -ve vlaue: Scan positive semi-axis and check there are no chars other than filling.
			for(int i= y0; i < strDepth.size()-1; i++){
				String x= strDepth.get(i);
				if(i == y0 && !x.equals(this.strFor1stNegU)){ // Only allowed char for zero line
					String msg= "Unexpected char in column bar. Got \"" 
							+ x + "\"; expected \"" + this.strFor1stNegU + "\". Bar:\n" + strDepth;
					throw new RuntimeException(msg);
				} else if(i > y0 && !x.equals(this.strForFill)){
					String msg= "Unexpected char in column bar. Got \"" + x 
							+ "\"; expected \"" + this.strForFill 
							+ "\". Bar:\n" + strDepth;
					throw new RuntimeException(msg);
				}
			}
		} else if (yValue >= 0) { // Score is positive or zero: Check negative semi-axis has only blanks 
			for(int i= 0; i <= y0; i++){
				String x= strDepth.get(i);
				if(i == y0){
					if(!x.equals(this.strForZero) && !x.equals(this.strFor1u) && !x.equals(this.strFor2u)){
						String msg= "Unexpected char in column bar. Got \"" + x + "\"" + "; Bar:\n" + strDepth;
						throw new RuntimeException(msg);
					}
				} else if(!x.equals(this.strForFill)){
					String msg= "Unexpected char in column bar. Got \"" + x + "\"" + "; Bar:\n" + strDepth;
					throw new RuntimeException(msg);
				}
				
			}
		}
	}

	
	@Deprecated
	public void TextProfileOld(List<Double> yValues, int yMaxLines){
		
		this.maxDepth= Integer.MIN_VALUE;;
        double minDepth= Integer.MAX_VALUE; 
		for(double y : yValues){
			if(y > this.maxDepth){
				this.maxDepth= y;
			}
			if(y < minDepth){
				minDepth= y;
			}
		}
		if(minDepth > 0){
			minDepth= 0;
		}

		this.scorePerDot= (double)(this.maxDepth - minDepth)/ (yMaxLines * 2); // * 2 because we use ':' for 2 units in a single line.
		// Locate zero on y axis:
		int y0= (int) Math.rint((double)(0 - minDepth) / (this.scorePerDot*2));
		List<List<String>> profile= new ArrayList<List<String>>();
		for(int i= 0; i < yValues.size(); i++){
			ArrayList<String> strDepth= new ArrayList<String>(); // This will be a vertical bar
			Double locDepth= (yValues.get(i) - minDepth) / this.scorePerDot; // yRescaled.get(i);
			
			if((int)Math.rint(locDepth) == 0){ // For zero coverage
				strDepth.add(this.strForZero); 
			}
			int nDouble= (int) ((Math.rint(locDepth)) / 2); // how many :
			for(int j= 0; j < nDouble; j++){
				strDepth.add(this.strFor2u);
			} 
			int diff= (int)Math.rint(locDepth- (nDouble*2)); // Get remainder
			if(diff == 2){
				strDepth.add(strFor2u);
			} else if(diff == 1) {
				strDepth.add(strFor1u);
			} else if(diff == 0){
				//
			} else {
				throw new RuntimeException("Unexpected division");
			}
			// Fill up list with blanks
			while(strDepth.size() < yMaxLines * 2){
				strDepth.add(strForFill);
		    }
			profile.add(strDepth);
		}
		this.profile= Utils.transpose(profile);
	}
	
	/*  G e t t e r s  */
	protected double getMaxDepth() {
		return maxDepth;
	}

	protected double getScorePerDot() {
		return scorePerDot;
	}

	protected List<List<String>> getProfile() {
		return profile;
	}

	protected String getStrForFill() {
		return strForFill;
	}

}
