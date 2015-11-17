package tracks;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.apache.commons.lang3.StringUtils;

import com.google.common.base.Joiner;

import htsjdk.samtools.BAMIndex;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import samTextViewer.GenomicCoords;
import samTextViewer.SamLocusIterator;
import samTextViewer.Utils;

@SuppressWarnings("deprecation")
public class TrackCoverage extends Track {

	/* A t t r i b u t e s */
	
	/** Max number of loci per window (column of text char). If you query 
	 * large intervals, only some loci will be stored. This const define the density of
	 * loci per window. */
	// final static int LOC_PER_WINDOW= 200;
	
	
	private List<ScreenLocusInfo> screenLocusInfoList= new ArrayList<ScreenLocusInfo>(); 

	/** Each dot in the screen output track corresponds to this many units of 
	 * score in the input. Typically this "reads per dot". */
	private double scorePerDot;   
	private double maxDepth;
	// private String inputBam;
	
	/* C o n s t r u c t o r */
	
	/**
	 * Construct coverage track from bam alignment in the provided interval. 
	 * Loci will be sampled according to the size of the interval and the size of the printable screen. 
	 * @param bam Input bam file
	 * @param gc Interval to sample positions from
	 * @param windowSize The size of the screen in number of characters.
	 * @param filters Record filters to apply to input sam records.
	 * @param bs Should loci be parsed also as BS-Seq data? 
	 * @throws IOException 
	 */
	public TrackCoverage(String bam, GenomicCoords gc,
			List<SamRecordFilter> filters, boolean bs) throws IOException{
		
		this.setGc(gc);
		this.setFilename(bam);
		SamReaderFactory srf=SamReaderFactory.make();
		srf.validationStringency(ValidationStringency.SILENT);
		SamReader samReader= srf.open(new File(bam));
		
		// SAMFileHeader fh= samReader.getFileHeader();
		//IntervalList il= locusSampler(fh, gc);
		IntervalList il= new IntervalList(samReader.getFileHeader());
		il.add(new Interval(gc.getChrom(), gc.getFrom(), gc.getTo()));
		SamLocusIterator samLocIter= new SamLocusIterator(samReader, il, true);
		samLocIter.setSamFilters(filters);
		Iterator<samTextViewer.SamLocusIterator.LocusInfo> iter= samLocIter.iterator();
	
		for(int i= 0; i < gc.getMapping().size(); i++){
			screenLocusInfoList.add(new ScreenLocusInfo());	
		}
		
		while(iter.hasNext()){			
			samTextViewer.SamLocusIterator.LocusInfo locusInfo= iter.next();
			int screenPos= Utils.getIndexOfclosestValue(locusInfo.getPosition(), gc.getMapping());
			byte refBase= '\0';
			if(gc.getRefSeq() != null){
				refBase= gc.getRefSeq()[screenPos];
			}
			screenLocusInfoList.get(screenPos).increment(locusInfo, refBase, bs);
		}
		samLocIter.close();
		samReader.close();
	}
	
	/* M e t h o d s */

	/**
	 * Printable coverage track. The height of the track in lines is `yMaxLines`.
	 * @param screenToGenomeMap List of genomic positions corresponding to each column on screen.
	 * @param yMaxLines
	 * @param rpm Should read counts be normalized by library size as Reads Per Million
	 * @return HashMapwith with keys/values the printable characteristics of the track. 
	 */
	@Override
	public String printToScreen(){
		
		if(this.getyMaxLines() == 0){return "";}
		
		List<Double> yValues= new ArrayList<Double>();
		for(ScreenLocusInfo x : screenLocusInfoList){
			yValues.add(x.getMeanDepth());
		}
		TextProfile textProfile= new TextProfile(yValues, this.getyMaxLines(), this.getYmin(), this.getYmax());
				
		this.maxDepth= textProfile.getMaxDepth();
		this.scorePerDot= textProfile.getScorePerDot();
		if(this.isRpm()){
			long libSize= getAlignedReadCount(new File(this.getFilename()));
			this.maxDepth= this.maxDepth / libSize * 1000000;
			this.scorePerDot= this.scorePerDot / libSize * 1000000;
		}
		
		ArrayList<String> lineStrings= new ArrayList<String>();
		for(int i= (textProfile.getProfile().size() - 1); i >= 0; i--){
			List<String> xl= textProfile.getProfile().get(i);
			Set<String> unique= new HashSet<String>(xl);
			lineStrings.add(StringUtils.join(xl, ""));
		}
		return Joiner.on("\n").join(lineStrings);
	}
        
    private long getAlignedReadCount(File bam){

    	SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);
		@SuppressWarnings("resource")
		BAMIndex sr= new SAMFileReader(bam).getIndex();
		long alnCount= 0; 
		int i= 0;
		while(true){
			try{
				alnCount += sr.getMetaData(i).getAlignedRecordCount();
			} catch(NullPointerException e){
				break;
			}
			i++;
		}
		sr.close();
		return alnCount;
    }
    
    /* S e t t e r s   and   G e t t e r s */
    
    /* This method makes sense calling only after having set the profile. Typically after  
     * */
    public double getScorePerDot(){
    	return this.scorePerDot;
    }
    public double getMaxDepth(){
    	this.printToScreen(); //  It's silly to call this just to set maxDepth.
    	return this.maxDepth;
    }    
    public List<ScreenLocusInfo> getScreenLocusInfoList(){
    	return screenLocusInfoList;
    }

	@Override
	public String getTitle(){
		return this.getFilename() + "; ylim: " + this.getYmin() + ", " + this.getYmax() + "; max: " + 
				Math.rint((this.getMaxDepth())*100)/100 + "; .= " + Math.rint((this.scorePerDot)) + ";\n";
	}
}
