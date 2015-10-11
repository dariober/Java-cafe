package samTextViewer;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;

/**
 * Ruler maps genome coordinates to screen coordinates. Be really sure you convert
 * genome intervals to screen usign the same algorithm though this project.
 * Coordinates are 1-based.
 * @author berald01
 */
public class Ruler {

	private int genomeFrom;
	private int genomeTo;
	private int screenWidth; 
	private List<Float> mapping= new ArrayList<Float>();
	
	/**
	 * Make sure Utils.compressNumericList() is consistent with splitting here.
	 * windowSize is the desired window size, this might differ from the effective 
	 * screenWidth.
	 */
	public Ruler(int genomeFrom, int genomeTo, int windowSize) {
		
		// First calculate the size of the screen window given the desired windowSize and 
		// genomic span from->to. The effective width on screen might differ from the desired
		// one because of the splitting. 
		
		// Once you have the effective screen width, just assign to each screen column the
		// corresponding genomic coordinate, rounded to nearest int.
		
		if(genomeFrom < 1 || genomeFrom > genomeTo){
			System.err.println("Invalid genome coordinates: from " + genomeFrom + " to " + genomeTo);
			System.exit(1);
		}
		
		this.genomeFrom= genomeFrom;
		this.genomeTo= genomeTo;
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
		}
		this.screenWidth= screenWidth;
		this.getMapping();
	}
	
	/* Methods */
	public float bpPerScreenColumn(){
		float bpPerScreenColumn= (genomeTo - genomeFrom + 1) / (float)this.screenWidth;
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
		int closest= Utils.getIndexOfclosestValue((float)genomePos, this.mapping);
/*		float bestDiff= Integer.MAX_VALUE;
		int closest= -1;
		for(int i= 0; i < this.mapping.size(); i++){ 
			// Iterate through entire list to find closest position on screen, it's a bit wasteful since
			// the list is ordered, but it's ok.
			float screenMark= this.mapping.get(i);
			float diff= Math.abs(genomePos - screenMark);
			if(diff < bestDiff){
				closest= i;
				bestDiff= diff;
			}
		}
		if(closest < 0){
			System.err.println("Invalid index position.");
			System.exit(1);
		} */
		return closest;
	}
	
	/**
	 * Mapping of genome position to each screen position.
	 */
	public  List<Float>  getMapping(){
		this.mapping.clear(); // This is important because if you run getMapping() multiple times you keep adding!
		for(int i= 0; i < this.screenWidth; i++){
			this.mapping.add(((i) * this.bpPerScreenColumn()) + this.genomeFrom);
		}
		return mapping;
	}
	
	public String toString(){
		String str= "Genome coords: " + this.genomeFrom + "-" + this.genomeTo 
				+ "; screen width: " + this.screenWidth
				+ "; scale: " + this.bpPerScreenColumn() + " bp/column" 
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
	
	public int getGenomeFrom() {
		return genomeFrom;
	}

	public int getGenomeTo() {
		return genomeTo;
	}

	public int getScreenWidth() {
		return screenWidth;
	}
}
