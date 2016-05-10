package tracks;

import java.util.ArrayList;
import java.util.List;

import samTextViewer.Utils;

/** Text representation of a continuous profile along the screen positions */
class TextProfile {
	
	private Double yMaxLimit= Double.NaN;   // Store the max depth of the track
	private Double yMinLimit= Double.NaN;
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
	 * @param yMimUser, yMaxUser Min and max values for y-axis as set by user. 
	 * If NaN use min and max from yValues  
	 */
	public TextProfile(List<Double> yValues, int yMaxLines, Double yMinUser, Double yMaxUser){
		
		// * Get ymin and ymax of input yValues
		double ymin= Double.MAX_VALUE;
		double ymax= Double.MIN_VALUE;
		for(double x : yValues){
			if(x > ymax){
				ymax= x;
			} 
			if(x < ymin){
				ymin= x;
			}
		}
		// * If yMinUser is NaN set it to ymin. Same for yMaxUser  
		if(yMinUser.isNaN()){
			yMinUser= ymin;
		}
		if(yMaxUser.isNaN()){
			yMaxUser= ymax;
		}
		if(yMinUser > 0){ // FIXME: min and max should be allowed to be <0
			yMinUser= 0.0;
		}
		if(yMaxUser < 0){
			yMaxUser= 0.0;
		}
		this.yMinLimit= yMinUser;
		this.yMaxLimit= yMaxUser;
		/*
		if(ymin.isNaN() != ymax.isNaN()){
			throw new RuntimeException("ymin and ymax must be both defined or both NaN");
		}
		if(ymin > ymax){
			throw new RuntimeException("Cannot have ymin >= ymax. Got ymin= " + ymin + "; ymax= " + ymax);
		}
		
		if(ymax.isNaN()){
			this.maxDepth= (double) Integer.MIN_VALUE;;
	        this.minDepth= (double) Integer.MAX_VALUE; 
			for(double y : yValues){
				if(y > this.maxDepth){
					this.maxDepth= y;
				}
				if(y < this.minDepth){
					this.minDepth= y;
				}
			}
			if(this.minDepth > 0){
				this.minDepth= 0.0;
			}
			if(this.maxDepth < 0){
				this.maxDepth= 0.0;
			}
		} else {
			this.maxDepth= ymax;
			this.minDepth= ymin;
		}
		*/
		this.scorePerDot= (double)(yMaxUser - yMinUser) / ((double)yMaxLines * 2); // * 2 because we use ':' for 2 units in a single line.
		
		// Locate zero on y axis. It's silly to generate a sequence just to find the index closest to zero. But anyway...
		List<Double> yAxis = Utils.seqFromToLenOut(yMinUser, yMaxUser, yMaxLines);
		int y0= Utils.getIndexOfclosestValue(0, yAxis);
		List<List<String>> profile= new ArrayList<List<String>>();
		for(int i= 0; i < yValues.size(); i++){
			double y= yValues.get(i);
			List<String> strDepth = prepareYColumn(y, yMaxLines, y0);
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
	private List<String> prepareYColumn(double yValue, int yMaxLines, int y0){
	
		Double yPosDotU= Math.abs(yValue / this.scorePerDot); // Y positions in dot units, not line units. 
		if((int)Math.rint(yPosDotU) == 0){ // For zero coverage
			ArrayList<String> strDepth= new ArrayList<String>();
			for(int j= 0; j < yMaxLines; j++){
				strDepth.add(strForFill);
			}			
			if(y0 < yMaxLines-1 || yMaxLines == 1){
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
			if(pos >= strDepth.size()) break;
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
		while(true){
			if(pos < 0) break;
			
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
				if(i == y0 && !(x.equals(this.strFor1stNegU) ||  x.equals(this.strForZero))){ // Only allowed char for zero line
					String msg= "y0= " + y0 + "; yValue= " + yValue +"; Unexpected char in column bar. Got \"" 
							+ x + "\"; expected \"" + this.strFor1stNegU + "\". Bar:\n" + strDepth;
					throw new RuntimeException(msg);
				} else if(i > y0 && !x.equals(this.strForFill)){
					String msg= "y0= " + y0 + "; yValue= " + yValue + "; Unexpected char in column bar. Got \"" + x 
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
						String msg= "y0= " + y0 + "; yValue= " + yValue + "; Unexpected char in column bar. Got \"" + x + "\"" + "; Bar:\n" + strDepth;
						throw new RuntimeException(msg);
					}
				} else if(!x.equals(this.strForFill)){
					String msg= "y0= " + y0 + "; yValue= " + yValue + "; Unexpected char in column bar. Got \"" + x + "\"" + "; Bar:\n" + strDepth;
					throw new RuntimeException(msg);
				}
				
			}
		}
	}
	
	/*  G e t t e r s  */
	protected double getYMaxLimit() {
		return yMaxLimit;
	}

	protected double getYMinLimit() {
		return yMinLimit;
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
