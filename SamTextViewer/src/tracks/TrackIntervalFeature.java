package tracks;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.lang3.StringUtils;

import samTextViewer.GenomicCoords;

public class TrackIntervalFeature {
 
	private List<IntervalFeature> intervalFeatureList= new ArrayList<IntervalFeature>();  
	private GenomicCoords gc; 
	
	/* Constructor */
	public TrackIntervalFeature(IntervalFeatureSet ifs, GenomicCoords gc) throws IOException {

		this.intervalFeatureList = ifs.getFeaturesInInterval(gc.getChrom(), gc.getFrom(), gc.getTo());
		for(IntervalFeature ift : intervalFeatureList){
			ift.mapToScreen(gc.getMapping());
		}
		this.gc= gc;
	}

	/* Methods */
	public String printToScreen(boolean noFormat) {
		
		String fwd= ">";
		String rev= "<";
		String na= "|";
		
		List<String> printable= new ArrayList<String>();
		for(int i= 0; i < gc.getMapping().size(); i++){ // First create empty track
			printable.add(" ");
		}
		for(IntervalFeature intervalFeature : intervalFeatureList){
			if(intervalFeature.getScreenFrom() == -1){
				continue; // Feature doesn't map to screen, this shouldn't happen though
			}
			for(int j= intervalFeature.getScreenFrom(); j <= intervalFeature.getScreenTo(); j++){
				String x= "";
				if(intervalFeature.getStrand() == '+'){
					x= (noFormat) ? fwd : "\033[48;5;147;38;5;240m" + fwd + "\033[0m";
				} else if(intervalFeature.getStrand() == '-'){
					x= (noFormat) ? rev : "\033[48;5;225;38;5;240m" + rev + "\033[0m";
				} else {
					x= (noFormat) ? na : "\033[48;5;250;38;5;240m" + na + "\033[0m";
				}
				printable.set(j, x);
			}
		}
		return StringUtils.join(printable, "");
	}
	
	public String toString(){
		StringBuilder sb= new StringBuilder();
		for(IntervalFeature ift : intervalFeatureList){
			sb.append(ift.getRaw());
			sb.append("\n");
		}
		return sb.toString().replaceAll("\n$", "");
	}

//	public String printGtfFeatures() {
//		StringBuilder sb= new StringBuilder();
//		for(IntervalFeature x : this.intervalFeatureList){
//			sb.append(x.toGtfString()).append("\n");
//		}
//		return sb.toString().replaceAll("\n$", "");
//	}	
}
