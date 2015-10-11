package coverageViewer;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPInputStream;

import org.apache.commons.lang3.StringUtils;

/**
 * Class containing sets of IntervalFeature objects. Essentially a representation of
 * a bed or gtf file. Similar to pybedtools BedTool object.
 * Essentially this is a HashMap with chrom as keys and List<IntervalFeature> as values.
 * Implementation is pretty naive but should be good enough.
 * @author berald01
 */
public class IntervalFeatureSet {
	
	private Map <String, List<IntervalFeature>> intervalMap= new HashMap <String, List<IntervalFeature>>(); 
	
	/* Constructor */
	/** 
	 * Construct from bed or gtf file.
	 * @throws IOException 
	 */
	public IntervalFeatureSet(File infile) throws IOException{
		
		Map <String, List<IntervalFeature>> intervalMap= new HashMap<String, List<IntervalFeature>>(); 
		
		BufferedReader br= null;
		InputStream gzipStream= null;
		if(infile.getName().endsWith(".gz")){
			InputStream fileStream = new FileInputStream(infile);
			gzipStream = new GZIPInputStream(fileStream);
			Reader decoder = new InputStreamReader(gzipStream, "UTF-8");
			br = new BufferedReader(decoder);
		} else {
			br = new BufferedReader(new FileReader(infile));		
		}
		String line;
		while ((line = br.readLine()) != null) {
			IntervalFeature f= new IntervalFeature(line, "bed");
			if(intervalMap.containsKey(f.getChrom())){
				intervalMap.get(f.getChrom()).add(f);
			} else {
				List<IntervalFeature> il= new ArrayList<IntervalFeature>();
				il.add(f);
				intervalMap.put(f.getChrom(), il);
			}
		}
		br.close();
		if(gzipStream != null){
			gzipStream.close();
		}
		this.intervalMap= intervalMap;
		this.sortIntervalsWithinChroms();
	}
	/* Methods */

	public List<IntervalFeature> getFeaturesInInterval(String chrom, int from, int to){
		if(from > to || from < 1 || to < 1){
			System.err.println("Invalid range: " + from + "-" + to);
			System.exit(1);
		}		
		
		List<IntervalFeature> xFeatures= new ArrayList<IntervalFeature>();

		List<IntervalFeature> thisChrom= this.getIntervalMap().get(chrom);		
		if(thisChrom == null){
			return xFeatures;
		}
		/* Overlap scenarios
		             from      to
		               |-------|
		     --------************-----  4.
		     ---------****------------- 1.
		     ------------****---------- 3.
		 	 -----------------****----- 2.
		     ----------------------***- Break iterating */
		for(IntervalFeature x : thisChrom){
			if( (x.getFrom() >= from && x.getFrom() <= to)  // 2. 3.
				|| (x.getTo() >= from && x.getTo() <= to)   // 1. 3.
				|| (x.getFrom() <= from && x.getTo() >= to)   // 4.
				){
				xFeatures.add(x);
			}
			if(x.getFrom() > to){
				break;
			}
		}
		return xFeatures;
	}
	
	private void sortIntervalsWithinChroms(){
		for(String chrom : this.intervalMap.keySet()){
			List<IntervalFeature> interalList = this.intervalMap.get(chrom);
			Collections.sort(interalList);
		}
	}
	
	/**
	 * Assign to all features in set the matching positions on screen. -1 if feature is
	 * not contained at all.
	 * */
	public void mapIntervalsToScreen(String chrom, List<Float> rulerMap){
		List<IntervalFeature> intervalsOnChrom= this.getIntervalMap().get(chrom);
		for(IntervalFeature f : intervalsOnChrom){
			f.mapToScreen(rulerMap);
		}
	}
	
	/**
	 * Return a string showing the mapping of this feature set to the given ruler.
	 * This is what the user sees.*/
	public String printableIntervals(String chrom, List<Float> rulerMap){
		// First prepare an empty list of length rulerMap. Each entry
		// Then map interval set to ruler and fill up the list where
		// a mapping is found.
		List<String> printable= new ArrayList<String>();
		for(int i= 0; i < rulerMap.size(); i++){
			printable.add("-");
		}
		this.mapIntervalsToScreen(chrom, rulerMap);
		List<IntervalFeature> mapSet = this.getIntervalMap().get(chrom);
		for(IntervalFeature x : mapSet){
			if(x.getScreenFrom() == -1){
				continue; // Feature doesn't map to screen
			}
			for(int j= x.getScreenFrom(); j == x.getScreenTo(); j++){
				if(x.getStrand() == '+'){
					printable.set(j, ">");
				} else if(x.getStrand() == '-'){
					printable.set(j, "<");
				} else {
					printable.set(j, "|");
				}
			}
		}
		return StringUtils.join(printable, "");
	}
	
	/* Setters and Getters */
	public Map<String, List<IntervalFeature>> getIntervalMap() {
		return intervalMap;
	}
	

	
}
