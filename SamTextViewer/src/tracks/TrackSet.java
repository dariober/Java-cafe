package tracks;

import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;

import org.apache.commons.lang3.text.StrTokenizer;

import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;

/** Class to hold tracks to be printed. 
 * */
public class TrackSet {
	
	private LinkedHashMap<String, Track> trackSet= new LinkedHashMap<String, Track>();
	
	/*   C o n s t r u c t o r s   */
	
	public TrackSet(){}
	
	/*   M e t h o d s   */

	public void addOrReplace(Track track){
		this.trackSet.put(track.getFileTag(), track);
	}
	
	/** From cmdInput extract regex and ylimits then iterate through the tracks list to set 
	 * the ylimits in the tracks whose filename matches the regex.
	 * The input list is updated in place! 
	*/
	public void setTrackYlimitsForRegex(String cmdInput) throws InvalidCommandLineException{

		StrTokenizer str= new StrTokenizer(cmdInput);
		str.setQuoteChar('\'');
		List<String> tokens= str.getTokenList();
		if(tokens.size() < 3){
			System.err.println("Error in ylim subcommand. Expected at least 3 args got: " + cmdInput);
			throw new InvalidCommandLineException();
		}
		String trackNameRegex= ".*"; // Default: Capture everything
		if(tokens.size() == 4){
			trackNameRegex= tokens.get(3);
		}
		
		try{
			Pattern.compile(trackNameRegex); // Validate regex
		} catch(PatternSyntaxException e){
	    	System.err.println("Invalid regex in: " + cmdInput);
	    	System.err.println(e.getDescription());
		}
		
		double ymin= Double.NaN;
        double ymax= Double.NaN;
		try{
			ymin= Double.parseDouble(tokens.get(1));
			ymax= Double.parseDouble(tokens.get(2));
		} catch(NumberFormatException e){
			ymin= Double.NaN;
			ymax= Double.NaN;
		}
		if(ymin >= ymax){
			System.err.println("Warning ymin >= ymax. Resetting to default.");
			ymin= Double.NaN;
			ymax= Double.NaN;							
		}
		for(Track tr : this.trackSet.values()){
			if(tr.getFileTag().matches(trackNameRegex)){
				tr.setYmin(ymin);
				tr.setYmax(ymax);
			}
		}
	}

	/** Set visibility for IntervalFeature tracks. 
	*/
	public void setVisibilityForTrackIntervalFeature(String cmdInput) throws InvalidCommandLineException{

		StrTokenizer str= new StrTokenizer(cmdInput);
		str.setQuoteChar('\'');
		List<String> tokens= str.getTokenList();

		// Defaults:
		String showRegex= ".*";  // Show all
		String hideRegex= "^$";    // Hide nothing
		String trackNameRegex= ".*"; // Apply to all tracks
		if(tokens.size() > 1){
			showRegex= tokens.get(1);
		}
		if(tokens.size() > 2){
			hideRegex= tokens.get(2);
		}
		if(tokens.size() > 3){
			trackNameRegex= tokens.get(3);
		}
		try{
			// Validate regex
			Pattern.compile(hideRegex); 
			Pattern.compile(showRegex);
			Pattern.compile(trackNameRegex); 
		} catch(PatternSyntaxException e){
	    	throw new PatternSyntaxException(e.getDescription(), cmdInput, -1);
		}

		System.err.println("Show: '" + showRegex + "'; hide: '" + hideRegex + "'; for tracks captured by '" + trackNameRegex + "':");
		for(Track tr : this.trackSet.values()){
			if(tr.getFileTag().matches(trackNameRegex)){
				System.err.println(tr.getFileTag());
				tr.setShowRegex(showRegex);
				tr.setHideRegex(hideRegex);
			}
		}
	}
	
	public GenomicCoords goToNextFeatureOnFile(String trackId, GenomicCoords curGc) throws InvalidGenomicCoordsException, IOException{

		trackId= trackId.trim();
		
		LinkedHashMap<String, Track> ifTracks = this.getIntervalFeatureTracks().getTrackSet();
		
		if(ifTracks.isEmpty()){
			System.err.println("\nTrackSet has no interval feature tracks!\n");
			return curGc;
		}
		
		Track tr;
		if(trackId.isEmpty() && ifTracks.size() == 1){
			tr= ifTracks.values().iterator().next();
		} else if(!ifTracks.containsKey(trackId)){ 	
			System.err.println("\nTag '" + trackId + "' not found in searchable track set:");
			System.err.println(ifTracks.keySet() + "\n");
			return curGc;
		} else {
			tr= ifTracks.get(trackId);
		}
		
		TrackIntervalFeature tif= (TrackIntervalFeature) tr;
		return tif.getIntervalFeatureSet().coordsOfNextFeature(curGc);
	}

	public GenomicCoords findNextStringOnTrack(String regex, String trackId, GenomicCoords curGc) throws InvalidGenomicCoordsException, IOException{

		trackId= trackId.trim();
		LinkedHashMap<String, Track> ifTracks = this.getIntervalFeatureTracks().getTrackSet();
		
		//if(ifTracks.isEmpty()){
		//	System.err.println("\nTrackSet has no interval feature tracks!\n");
		//	return curGc;
		//}
		
		Track tr;
		if(trackId.isEmpty() && ifTracks.size() == 1){
			tr= ifTracks.values().iterator().next();
		} else if(!ifTracks.containsKey(trackId)){ 	
			System.err.println("\nTag '" + trackId + "' not found in searchable track set:");
			System.err.println(ifTracks.keySet() + "\n");
			return curGc;
		} else {
			tr= ifTracks.get(trackId);
		}
		
		TrackIntervalFeature tif= (TrackIntervalFeature) tr;
		return tif.getIntervalFeatureSet().findNextRegex(curGc, regex);
	}

	private TrackSet getIntervalFeatureTracks(){
		TrackSet ifSet= new TrackSet();
		for(Track tr : this.trackSet.values()){
			if(Utils.getFileTypeFromName(tr.getFilename()).equals(TrackFormat.BED) 
			   || Utils.getFileTypeFromName(tr.getFilename()).equals(TrackFormat.GFF)){
				ifSet.addOrReplace(tr);
			}
		}
		return ifSet;
	}
	
	public void selectDataColumnForBedgraph(int bdgDataColIdx, String trackIdRegex){
		
		for(Track tr : this.trackSet.values()){
			
			if(Utils.getFileTypeFromName(tr.getFilename()).equals(TrackFormat.BEDGRAPH) &&
			   tr.getFileTag().matches(trackIdRegex)){
		
				TrackWiggles bdg= (TrackWiggles) tr;
				bdg.setBdgDataColIdx(bdgDataColIdx);
			
			}
		}		
	}
	
	// STUB:
	//public void setTrackHeight(int height, String trackIdRegex){
	//	for(Track tr : this.trackSet.values()){
	//		if(tr.getFileTag().matches(trackIdRegex)){
	//			tr.setyMaxLines(height);
	//		}
	//	}		
	//}
	
	/*   S e t t e r s   and   G e t t e r s  */
	public LinkedHashMap<String, Track> getTrackSet() {
		return trackSet;
	}
	
}
