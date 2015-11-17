package tracks;

import java.util.ArrayList;
import java.util.List;

import samTextViewer.GenomicCoords;

public class Track {

	private String title= "";
	private int yMaxLines= 10;
	private String filename= "N/A";
	private List<Double> screenScores= new ArrayList<Double>();
	private GenomicCoords gc;
	private boolean rpm= false;
	private boolean noFormat= false; 
	private double ymin= Double.NaN; // Same as R ylim()
	private double ymax= Double.NaN;

//	public Track(){}
	
	/* Printers */
	public String printToScreen(){
		return null;
	}

	public String printFeatures(){
		return "";
	};
	
	public String toString(){
		return this.getFilename();
	}
	
	/* Setters and getters */
	public void setTitle(String title){
		this.title= title;
	}
	public String getTitle(){
		return this.title;
	}
	public int getyMaxLines() {
		return yMaxLines;
	}
	public void setyMaxLines(int yMaxLines) {
		this.yMaxLines = yMaxLines;
	}
	public String getFilename() {
		return filename;
	}
	public void setFilename(String filename) {
		this.filename = filename;
	}
	protected List<Double> getScreenScores() {
		return screenScores;
	}
	protected void setScreenScores(List<Double> screenScores) {
		this.screenScores = screenScores;
	}

	protected GenomicCoords getGc() {
		return gc;
	}

	protected void setGc(GenomicCoords gc) {
		this.gc = gc;
	}

	public boolean isRpm() {
		return rpm;
	}

	public void setRpm(boolean rpm) {
		this.rpm = rpm;
	}

	public boolean isNoFormat() { return noFormat; }
	public void setNoFormat(boolean noFormat) { this.noFormat = noFormat; }

	public double getYmin() { return ymin;}
	public void setYmin(double ymin) { this.ymin = ymin;}

	public double getYmax() { return ymax; }
	public void setYmax(double ymax) { this.ymax = ymax; }

}

