package tracks;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.lang3.StringUtils;

import samTextViewer.GenomicCoords;

public class TrackIntervalFeature extends Track {
 
	private List<IntervalFeature> intervalFeatureList= new ArrayList<IntervalFeature>();  
	private IntervalFeatureSet intervalFeatureSet;
	
	/* Constructor */
	
	public TrackIntervalFeature(String filename, GenomicCoords gc) throws IOException{
		this.setGc(gc);
		this.setFilename(filename);
		this.intervalFeatureSet= new IntervalFeatureSet(new File(filename));
		this.update();
	}
	
	/* Methods */
	
	public void update() throws IOException{
		this.intervalFeatureList = this.intervalFeatureSet.getFeaturesInInterval(
				this.getGc().getChrom(), this.getGc().getFrom(), this.getGc().getTo());
		for(IntervalFeature ift : intervalFeatureList){
			ift.mapToScreen(this.getGc().getMapping());
		}
	}
	
	@Override
	public String printToScreen() {
		
		String fwd= ">";
		String rev= "<";
		String na= "|";
		
		List<String> printable= new ArrayList<String>();
		for(int i= 0; i < this.getGc().getMapping().size(); i++){ // First create empty track
			printable.add(" ");
		}
		for(IntervalFeature intervalFeature : intervalFeatureList){
			if(intervalFeature.getScreenFrom() == -1){
				continue; // Feature doesn't map to screen, this shouldn't happen though
			}
			for(int j= intervalFeature.getScreenFrom(); j <= intervalFeature.getScreenTo(); j++){
				String x= "";
				if(intervalFeature.getStrand() == '+'){
					x= (this.isNoFormat()) ? fwd : "\033[48;5;147;38;5;240m" + fwd + "\033[0m";
				} else if(intervalFeature.getStrand() == '-'){
					x= (this.isNoFormat()) ? rev : "\033[48;5;225;38;5;240m" + rev + "\033[0m";
				} else {
					x= (this.isNoFormat()) ? na : "\033[48;5;250;38;5;240m" + na + "\033[0m";
				}
				printable.set(j, x);
			}
		}
		return StringUtils.join(printable, "");
	}
	
	@Override
	public String printFeatures(){
		StringBuilder sb= new StringBuilder();
		for(IntervalFeature ift : intervalFeatureList){
			sb.append(ift.getRaw());
			sb.append("\n");
		}
		return sb.toString(); // NB: Leave last trailing /n
	}
	
	@Override
	public String getTitle(){
		return this.getFileTag() + ";\n";
	}
	
	protected IntervalFeatureSet getIntervalFeatureSet() { return intervalFeatureSet; }

}
