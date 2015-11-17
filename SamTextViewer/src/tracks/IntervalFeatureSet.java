package tracks;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
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
import java.util.Set;
import java.util.zip.GZIPInputStream;

import exceptions.InvalidGenomicCoordsException;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.tribble.readers.TabixReader.Iterator;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;

/**
 * Class containing sets of IntervalFeature objects. Essentially a representation of
 * a bed or gtf file. Similar to pybedtools BedTool object.
 * Essentially this is a HashMap with chrom as keys and List<IntervalFeature> as values.
 * Implementation is pretty naive but should be good enough.
 * @author berald01
 */
public class IntervalFeatureSet {
	
	private Map <String, List<IntervalFeature>> intervalMap; // new HashMap <String, List<IntervalFeature>>(); 
	private TabixReader tabixReader= null;
	private boolean isTabix= false;
	private String type;
	
	/* Constructor */
	/** Construct from bed or gtf file.
	 * @throws IOException */
	public IntervalFeatureSet(File infile) throws IOException{
		
		this.type= Utils.getFileTypeFromName(infile.getName());
		
		if(new File(infile.getAbsolutePath() + ".tbi").exists()){
			this.tabixReader= new TabixReader(infile.getAbsolutePath());
			this.isTabix= true;
		} else {
			this.intervalMap= loadFileIntoIntervalMap(infile);
			this.sortIntervalsWithinChroms();
			this.isTabix= false;
		}
	}
	/* Methods */

	public List<IntervalFeature> getFeaturesInInterval(String chrom, int from, int to) throws IOException{
		if(from > to || from < 1 || to < 1){
			throw new RuntimeException("Invalid range: " + from + "-" + to);
		}		
		List<IntervalFeature> xFeatures= new ArrayList<IntervalFeature>();
		if(isTabix){
			Iterator qry = this.tabixReader.query(chrom,  from, to);
			while(true){
				String q = qry.next();
				if(q == null){
					break;
				}
				IntervalFeature intervalFeature= new IntervalFeature(q, this.type);
				xFeatures.add(intervalFeature);
			}
			return xFeatures;
		} else {
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
	}
	
	private Map <String, List<IntervalFeature>> loadFileIntoIntervalMap(File infile) throws IOException{
		
		System.err.print("Reading file '" + infile.getName() + "'...");
		
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
			if(line.trim().startsWith("#")){
				continue;
			}
			IntervalFeature f= new IntervalFeature(line, this.type);
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
		System.err.println(" Done");
		return intervalMap;
	} 
	
	private void sortIntervalsWithinChroms(){
		for(String chrom : this.intervalMap.keySet()){
			List<IntervalFeature> interalList = this.intervalMap.get(chrom);
			Collections.sort(interalList);
		}
	}
	
	/* Setters and Getters */
	public Map<String, List<IntervalFeature>> getIntervalMap() {
		return intervalMap;
	}

	/** Get the next feature on chrom after "from" position or null if no 
	 * feature found 
	 * @throws IOException */
	private IntervalFeature getNextFeatureOnChrom(String chrom, int from) throws IOException{
		
		if(this.intervalMap != null){
			List<IntervalFeature> featuresList = this.intervalMap.get(chrom);
			if(featuresList == null){ // No features at all on this chrom.
				return null;
			}
			for(IntervalFeature x : featuresList){
				if(x.getFrom() > from){
					return x;
				}
			}
			return null;
		} else if(this.isTabix){
			Iterator iter = this.tabixReader.query(chrom, from, Integer.MAX_VALUE);
			while(true){
				String line= iter.next();
				if(line == null){
					return null;
				} 
				IntervalFeature x= new IntervalFeature(line, this.type);
				if(x.getFrom() > from){
					return x;
				}
			}
		}
		return null;
	}
	
	/** Return the set chroms sorted but with and first chrom set to startChrom.
	 *  */
	protected List<String> getChromListStartingAt(Set<String> chroms, String startChrom){
		// Set:            chr1 chr2 chr3 chr4 chr5 chr6 chr7
		// StartChrom:     chr3
		// Ordered chroms: chr3 chr4 chr5 chr6 chr7 chr1 chr2
		List<String> orderedChroms= new ArrayList<String>();
		orderedChroms.addAll(chroms);
		Collections.sort(orderedChroms);
		int idx= orderedChroms.indexOf(startChrom);
		if(idx == -1){ // If startChrom is not present at all in the bed/gtf file.
			return orderedChroms;
		}
		List<String> chromsStartingAt= new ArrayList<String>();
		chromsStartingAt.addAll(orderedChroms.subList(idx, orderedChroms.size()));
		chromsStartingAt.addAll(orderedChroms.subList(0, idx));
		return chromsStartingAt;
	}
	
	/** Search chrom to find the *next* feature matching the given string. The search will 
	 * wrap around the chrom if not found in the chunk following "from". */
	protected IntervalFeature findNextStringOnChrom(String string, String chrom, int from) throws IOException{
		
		if(this.intervalMap != null){
	
			// We start search from input chrom and starting position. We'll return to 
			// the start of this chrom only after having searched all the other chroms.
			int startingPoint= from;
			List<String> chromSearchOrder = getChromListStartingAt(this.intervalMap.keySet(), chrom);
			chromSearchOrder.add(chrom);
			for(String curChrom : chromSearchOrder){
				
				List<IntervalFeature> featuresList = this.intervalMap.get(curChrom);
				for(IntervalFeature x : featuresList){
					if(x.getFrom() > startingPoint && x.getRaw().toLowerCase().contains(string.toLowerCase())){
						return x;
					}
				}
				startingPoint= 0;
				
			} return null; // Not found anywhere

		} else if(this.isTabix){
			
			int startingPoint= from;
			List<String> chromSearchOrder = getChromListStartingAt(this.tabixReader.getChromosomes(), chrom);
			chromSearchOrder.add(chrom);
			for(String curChrom : chromSearchOrder){
				
				Iterator iter = this.tabixReader.query(curChrom , startingPoint, Integer.MAX_VALUE);
				while(true){
					String line= iter.next();
					if(line == null) break;
					if(line.toLowerCase().contains(string.toLowerCase())){
						IntervalFeature x= new IntervalFeature(line, this.type);
						if(x.getFrom() > startingPoint){
							return x;
						}
					}
				}
				startingPoint= 0;
			} return null; // Not found anywhere
		}
		return null;
	}
	
	public GenomicCoords findNextString(GenomicCoords currentGc, String string) throws IOException, InvalidGenomicCoordsException{

		IntervalFeature nextFeature= findNextStringOnChrom(string, currentGc.getChrom(), currentGc.getTo());
		if(nextFeature == null){
			return currentGc;
		}
		GenomicCoords nextGc= new GenomicCoords(
				nextFeature.getChrom(), 
				nextFeature.getFrom(), 
				nextFeature.getFrom() + currentGc.getGenomicWindowSize() - 1, 
				currentGc.getSamSeqDict(),
				currentGc.getUserWindowSize(),
				currentGc.getFastaFile());
		return nextGc;
	}
	
	public GenomicCoords coordsOfNextFeature(GenomicCoords currentGc) throws InvalidGenomicCoordsException, IOException {
		IntervalFeature nextFeature= getNextFeatureOnChrom(currentGc.getChrom(), currentGc.getTo());
		if(nextFeature == null){
			return currentGc;
		}
		GenomicCoords nextGc= new GenomicCoords(
				nextFeature.getChrom(), 
				nextFeature.getFrom(), 
				nextFeature.getFrom() + currentGc.getGenomicWindowSize() -1, 
				currentGc.getSamSeqDict(),
				currentGc.getUserWindowSize(),
				currentGc.getFastaFile());
		return nextGc;
	}
}
