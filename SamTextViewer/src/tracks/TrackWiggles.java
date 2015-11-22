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
import java.util.List;
import java.util.zip.GZIPInputStream;

import org.apache.commons.lang3.StringUtils;
import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;
import org.broad.igv.tdf.TDFUtils;

import com.google.common.base.Joiner;

import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndex;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.tribble.readers.TabixReader.Iterator;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;

/** Process wiggle file formats. Mostly using IGV classes. 
 * bigBed, bigWig, */
public class TrackWiggles extends Track {

	private double maxDepth;
	private double scorePerDot;
	private List<ScreenWiggleLocusInfo> screenWiggleLocusInfoList;
	private int bdgDataColIdx= 4; 
	
	/* C o n s t r u c t o r s */

	/**
	 * Read bigWig from local file or remote URL.
	 * @param filename Filename or URL to access 
	 * @param gc Query coordinates and size of printable window 
	 * @throws IOException */
	public TrackWiggles(String filename, GenomicCoords gc, int bdgDataColIdx) throws IOException{

		this.setGc(gc);
		this.setFilename(filename);
		this.bdgDataColIdx= bdgDataColIdx;
		this.update();
		
	};
	
	public void update() throws IOException {

		if(this.bdgDataColIdx < 4){
			System.err.println("Invalid index for bedgraph column of data value. Resetting to 4. Expected >=4. Got " + this.bdgDataColIdx);
			this.bdgDataColIdx= 4;
		}

		if(Utils.getFileTypeFromName(this.getFilename()).equals(TrackFormat.BIGWIG)){
			
			BBFileReader reader=new BBFileReader(this.getFilename()); // or url for remote access.
			if(!reader.getBBFileHeader().isBigWig()){
				System.err.println("Invalid file type " + this.getFilename());
				System.exit(1);			
			}
			bigWigToScores(reader);
			
		} else if(Utils.getFileTypeFromName(this.getFilename()).equals(TrackFormat.TDF)){

			this.screenWiggleLocusInfoList= 
					TDFUtils.tdfRangeToScreen(this.getFilename(), this.getGc().getChrom(), 
							this.getGc().getFrom(), this.getGc().getTo(), this.getGc().getMapping());
			
			ArrayList<Double> screenScores= new ArrayList<Double>();
			for(ScreenWiggleLocusInfo x : screenWiggleLocusInfoList){
				screenScores.add((double)x.getMeanScore());
			}
			this.setScreenScores(screenScores);	
			
		} else if(Utils.getFileTypeFromName(this.getFilename()).equals(TrackFormat.BEDGRAPH)){

			if(Utils.hasTabixIndex(this.getFilename())){
				bedGraphToScores(this.getFilename());
			} else if(Utils.hasTabixIndex(this.getFilename() + ".samTextViewer.tmp.gz")){
				bedGraphToScores(this.getFilename() + ".samTextViewer.tmp.gz");
			} else {
				blockCompressAndIndex(this.getFilename(), this.getFilename() + ".samTextViewer.tmp.gz", true);
				bedGraphToScores(this.getFilename() + ".samTextViewer.tmp.gz");
			}
		} else {
			throw new RuntimeException("Extension (i.e. file type) not recognized for " + this.getFilename());
		}
	}

	
	/*  M e t h o d s  */
	/**
	 * Block compress input file and create associated tabix index. Newly created file and index are
	 * deleted on exit if deleteOnExit true.
	 * @throws IOException 
	 * */
	private void blockCompressAndIndex(String in, String bgzfOut, boolean deleteOnExit) throws IOException {
		
		BufferedReader br= null;
		InputStream gzipStream= null;
		if(in.endsWith(".gz")){
			InputStream fileStream = new FileInputStream(in);
			gzipStream = new GZIPInputStream(fileStream);
			Reader decoder = new InputStreamReader(gzipStream, "UTF-8");
			br = new BufferedReader(decoder);
		} else {
			br = new BufferedReader(new FileReader(in));
		}
		System.err.print("Block compressing " + in + "... ");
		BlockCompressedOutputStream blockOs= new BlockCompressedOutputStream(bgzfOut);
		String line;
		boolean isFirst= true;
		while ((line = br.readLine()) != null) {
			if(line.trim().startsWith("#")){
				continue;
			}
			if(isFirst){
				isFirst= false;
				if(!isValidBedGraphLine(line)){ // Allow first line to fail: Might be header.
					System.err.print("First line skipped. ");
					continue;
				}
			}
			line += "\n";
			blockOs.write(line.getBytes());
		}
		br.close();
		if(gzipStream != null){
			gzipStream.close();
		}
		blockOs.close();
		
		System.err.print("Indexing... ");
		File bgzfFile= new File(bgzfOut);
		BEDCodec codec=new BEDCodec();
		TabixIndex tabixIndex =
				IndexFactory.createTabixIndex(bgzfFile, codec, TabixFormat.BED, null);
		tabixIndex.writeBasedOnFeatureFile(bgzfFile);
		System.err.println("Done");
		if(deleteOnExit){
			bgzfFile.deleteOnExit();
			File idx= new File(bgzfFile + ".tbi");
			idx.deleteOnExit();
		}
	}

	/** Return true if line looks like a valid bedgraph record  
	 * */
	public static boolean isValidBedGraphLine(String line){
		String[] bdg= line.split("\t");
		if(bdg.length < 4){
			return false;
		}
		try{
			Integer.parseInt(bdg[1]);
			Integer.parseInt(bdg[2]);
		} catch(NumberFormatException e){
			return false;
		}
		return true;
	}
	
	/** Populate object using bigWig data */
	private void bigWigToScores(BBFileReader reader){

		// List of length equal to screen size. Each inner map contains info about the screen locus 
		List<ScreenWiggleLocusInfo> screenWiggleLocusInfoList= new ArrayList<ScreenWiggleLocusInfo>();
		for(int i= 0; i < getGc().getUserWindowSize(); i++){
			screenWiggleLocusInfoList.add(new ScreenWiggleLocusInfo());
		}
		
		BigWigIterator iter = reader.getBigWigIterator(getGc().getChrom(), getGc().getFrom(), getGc().getChrom(), getGc().getTo(), false);
		while(iter.hasNext()){
			WigItem bw = iter.next();
			for(int i= bw.getStartBase(); i <= bw.getEndBase(); i++){
				int idx= Utils.getIndexOfclosestValue(i, getGc().getMapping()); // Where should this position be mapped on screen?
				screenWiggleLocusInfoList.get(idx).increment(bw.getWigValue());
			} 
		}
		reader.close();
		ArrayList<Double> screenScores= new ArrayList<Double>();
		for(ScreenWiggleLocusInfo x : screenWiggleLocusInfoList){
			screenScores.add((double)x.getMeanScore());
		}
		this.setScreenScores(screenScores);		
	}
	
	/** Get values for bedgraph
	 * */
	private void bedGraphToScores(String fileName) throws IOException{
		
		List<ScreenWiggleLocusInfo> screenWiggleLocusInfoList= new ArrayList<ScreenWiggleLocusInfo>();
		for(int i= 0; i < getGc().getUserWindowSize(); i++){
			screenWiggleLocusInfoList.add(new ScreenWiggleLocusInfo());
		}
		
		try {
			TabixReader tabixReader= new TabixReader(fileName);
			Iterator qry= tabixReader.query(this.getGc().getChrom(), this.getGc().getFrom()-1, this.getGc().getTo());
			while(true){
				String q = qry.next();
				if(q == null){
					break;
				}
				String[] tokens= q.split("\t");
				int screenFrom= Utils.getIndexOfclosestValue(Integer.valueOf(tokens[1])+1, this.getGc().getMapping());
				int screenTo= Utils.getIndexOfclosestValue(Integer.valueOf(tokens[2]), this.getGc().getMapping());
				float value= Float.valueOf(tokens[this.bdgDataColIdx-1]);
				for(int i= screenFrom; i <= screenTo; i++){
					screenWiggleLocusInfoList.get(i).increment(value);
				}
			}
		} catch (IOException e) {			
			e.printStackTrace();
			System.err.println("Could not open tabix file: " + fileName);
			System.err.println("Is the file sorted and indexed? After sorting by position (sort e.g. -k1,1 -k2,2n), compress with bgzip and index with e.g.:");
			System.err.println("\nbgzip " + fileName);
			System.err.println("tabix -p bed " + fileName + "\n");
		}
		ArrayList<Double> screenScores= new ArrayList<Double>();
		for(ScreenWiggleLocusInfo x : screenWiggleLocusInfoList){
			screenScores.add((double)x.getMeanScore());
		}
		this.setScreenScores(screenScores);
		return;
	}
	
	@Override
	public String printToScreen(){ // int yMaxLines, Double ymin, Double ymax
		
		if(this.getyMaxLines() == 0){return "";}
		
		TextProfile textProfile= new TextProfile(this.getScreenScores(), this.getyMaxLines(), this.getYmin(), this.getYmax());
		
		this.scorePerDot= textProfile.getScorePerDot();
		this.maxDepth= textProfile.getMaxDepth();

		ArrayList<String> lineStrings= new ArrayList<String>();
		for(int i= (textProfile.getProfile().size() - 1); i >= 0; i--){
			List<String> xl= textProfile.getProfile().get(i);
			lineStrings.add(StringUtils.join(xl, ""));
		}
		return Joiner.on("\n").join(lineStrings);
	}
		
	/*   S e t t e r s   and   G e t t e r s */
	
	public double getMaxDepth() {
		return maxDepth;
	}

	public double getScorePerDot() {
		return scorePerDot;
	}

	protected int getBdgDataColIdx() { return bdgDataColIdx; }
	protected void setBdgDataColIdx(int bdgDataColIdx) { this.bdgDataColIdx = bdgDataColIdx; }

	@Override
	public String getTitle(){
		return this.getFileTag() + "; ylim: " + this.getYmin() + ", " + this.getYmax() + "; max: " + 
				Math.rint((this.maxDepth)*100)/100 + "; .= " + Math.rint((this.scorePerDot)*100)/100 + ";\n";
	}
	
}
