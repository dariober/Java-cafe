package samTextViewer;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;

/**
 * Ruler_TO_BE_DEPRECTED maps genome coordinates to screen coordinates. Be really sure you convert
 * genome intervals to screen usign the same algorithm though this project.
 * Coordinates are 1-based.
 * @author berald01
 */
public class Ruler_TO_BE_DEPRECTED {

	/** List extending from genomeFrom to genomeTo and of length `screenWidth` or less. Each
	 * index of this list is an index of the column on screen. This list is to be used to map features 
	 * (coverage, bed intervals etc) to the screen 
	 * */
	private List<Double> mapping= new ArrayList<Double>();

	private int genomeFrom; // Internal use only. Do not create setters and getters.
	private int genomeTo; // Internal use only. Do not create setters and getters.
	private int screenWidth; // Internal use only. Do not create setters and getters.
	

	/** 
	 * Mapping of genomic coordinates to screen coordinates.
	 * */
	public Ruler_TO_BE_DEPRECTED(int genomeFrom, int genomeTo, int windowSize) {
				
		if(genomeFrom < 1 || genomeFrom > genomeTo){
			System.err.println("Invalid genome coordinates: from " + genomeFrom + " to " + genomeTo);
			System.exit(1);
		}
		
		this.mapping= this.seqFromToLenOut(genomeFrom, genomeTo, windowSize);
		this.genomeFrom= genomeFrom;
		this.genomeTo= genomeTo;
		this.screenWidth= this.mapping.size();
	}
	
	/* Methods */
	
	/**
	 * Same as R seq(from to, length.out). 
	 * @param from
	 * @param to
	 * @param lengthOut Length of sequence, effectively the desired screen width.
	 * @return
	 */
	private List<Double> seqFromToLenOut(int from, int to, int lengthOut){
		
		if(from < 1 || from > to){
			System.err.println("Invalid genome coordinates: from " + from + " to " + to);
			System.exit(1);
		}
		this.genomeFrom= from;
		this.genomeTo= to;
	
		int span= to - from + 1;
		// If the genomic span is less then screen size, reduce screen size to.
		// If genomic span == screenSize then you have a mapping one to one.
		if(span <= lengthOut){ 
			for(int i= from; i <= to; i++){
				this.mapping.add((double)i);
			}
			return this.mapping;
		}
		
		double step= ((double)span - 1)/(lengthOut - 1);
		this.mapping.add((double)from);
		for(int i= 1; i < lengthOut; i++){
			this.mapping.add((double)this.mapping.get(i-1)+step);
		}
		
		// First check last point is close enough to expection. If so, replace last point with
		// exact desired.
		double diffTo= Math.abs(this.mapping.get(this.mapping.size() - 1) - to);
		if(diffTo > ((float)to * 0.001)){
			System.err.println("Error generating sequence:");
			System.err.println("Last point: " + this.mapping.get(this.mapping.size() - 1));
			System.err.println("To diff: " + diffTo);
			System.err.println("Step: " + step);
			System.err.println("Length: " + mapping.subList(33320, 33332));
		} else {
			this.mapping.set(mapping.size()-1, (double)to);
		}
		
		double diffFrom= Math.abs(this.mapping.get(0) - from);		
		if(diffFrom > 0.01 || this.mapping.size() != lengthOut){
			System.err.println("Error generating sequence:");
			System.err.println("Expected size: " + lengthOut + "; Effective: " + this.mapping.size());
			System.err.println("From diff: " + diffFrom);
			System.exit(1);
		}
		return this.mapping;
	}
	
	public float getBpPerScreenColumn(){
		float bpPerScreenColumn= (this.genomeTo - this.genomeFrom + 1) / screenWidth;
		return bpPerScreenColumn;
	}

	/**
	 * Mapping of genome positions to screen. Screen position is
	 * 0-based. I.e. first column has position 0. For genomic
	 * spans only slightly larger then screen size the mapping might be
	 * inaccurate due to rounding during preparation of windows. 
	 * */
	public int getScreenPositionAtGenomePosition(int genomePos){
		
		if(genomePos < this.genomeFrom || genomePos > this.genomeTo){
			return -1;
		}
		int closest= Utils.getIndexOfclosestValue((double)genomePos, this.mapping);
		return closest;
	}
	
	public String toString(){
		String str= "Genome coords: " + this.genomeFrom + "-" + this.genomeTo 
				+ "; screen width: " + this.screenWidth
				+ "; scale: " + this.getBpPerScreenColumn() + " bp/column" 
				+ "; Mapping: " + this.mapping;
		return str;
	}
	
	public String printableRuler(int markDist){
    	String numberLine= "";
    	int prevLen= 0;
    	int i= 0;
		while(i < this.screenWidth){
			String posMark= String.valueOf(Math.round(this.mapping.get(i)));
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

	/* Setters and Getters */
	public List<Double> getMapping() {
		return mapping;
	}	
}

/*
int genomeSpan= genomeTo - genomeFrom + 1;
int grpSize= Math.round(genomeSpan / windowSize);
int screenWidth= 0;
int sublistLen= 0; // Counterpart of sublist in Utils.compressNumericList()
for(int i= 0; i < genomeSpan; i++){
	sublistLen++;
	if(sublistLen == grpSize || grpSize < 1){ // grpSize < 1 means exact mapping of 1bp = text column 
		sublistLen= 0;
		screenWidth++;
	}
}
if(sublistLen > 0){
	screenWidth++;			
} */
// this.screenWidth= screenWidth;
// this.getMapping();

/**TO BE DEPRECATED
 * Mapping of genome position to each screen position.
 */
/*public  List<Float>  getMapping(){
	this.mapping.clear(); // This is important because if you run getMapping() multiple times you keep adding!
	for(int i= 0; i < this.screenWidth; i++){
		this.mapping.add(((i) * this.bpPerScreenColumn()) + this.genomeFrom);
	}
	return mapping;
} */

