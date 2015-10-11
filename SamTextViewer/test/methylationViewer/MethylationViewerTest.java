package methylationViewer;

import filter.FlagToFilter;
import filter.ReadFromTopStrandFilter;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.SamLocusIterator;
import htsjdk.samtools.util.SamLocusIterator.LocusInfo;
import htsjdk.samtools.util.SamLocusIterator.RecordAndOffset;

import java.util.List;

import org.junit.Test;

import readWriteBAMUtils.ReadWriteBAMUtils;

public class MethylationViewerTest {

	@Test
	public void canReturnMethylCountsForInterval(){

		SamReader samReader= ReadWriteBAMUtils.reader("test_data/mjb050_oxBS.bam", ValidationStringency.SILENT);
		SAMFileHeader fh= samReader.getFileHeader();
		IntervalList il= new IntervalList(fh);
		Interval interval= new Interval("chrY", 1, 100);
		il.add(interval);		
		List<SamRecordFilter> filters= FlagToFilter.flagToFilterList(0, 0);
		
		SamLocusIterator samLocIter= new SamLocusIterator(samReader, il, true);
		samLocIter.setSamFilters(filters);
		LocusInfo locus= samLocIter.next();
		
		
		
	}
	
	@Test
	public void STUBcanCallBS(){
		
		
		int f_incl= 0;
		int F_excl= 0;
		String chrom= "chrY";
		int from= 1;
		int to= 1;
		
		SamReader samReader= ReadWriteBAMUtils.reader("test_data/mjb050_oxBS.bam", ValidationStringency.SILENT);
		SAMFileHeader fh= samReader.getFileHeader();
		IntervalList il= new IntervalList(fh);
		Interval interval= new Interval(chrom, from, to);
		il.add(interval);
		
		List<SamRecordFilter> filters= FlagToFilter.flagToFilterList(f_incl, F_excl);
	
		SamLocusIterator samLocIter= new SamLocusIterator(samReader, il, true);
		samLocIter.setSamFilters(filters);
		while(samLocIter.hasNext()){
			LocusInfo locus= samLocIter.next();
			int M= 0;
			int U= 0;
			int mism= 0;
			for(RecordAndOffset recOff : locus.getRecordAndPositions()){
				int pos= locus.getPosition();
				// Code to get ref sequence at pos
				
				// If ref sequence is C, count read bases if:
				char refbase= Character.toUpperCase('\0');
				char rb= Character.toUpperCase((char)recOff.getReadBase());
				
				boolean isTopStrand= !(new ReadFromTopStrandFilter(true)).filterOut(recOff.getRecord());
				
				if(refbase == 'C'){

					if( isTopStrand	){ // -ve 2nd pair
						if(rb == 'C'){
							M++;
						} else if(rb == 'T'){
							U++;
						} else {
							mism++;
						}
					}  					
				} else if (refbase == 'G'){

					if(	!isTopStrand ){
						if(rb == 'G'){
							M++;
						} else if(rb == 'A'){
							U++;
						} else {
							mism++;
						}
							// System.out.println(locus.getPosition() + " ");					
						}  
				} else {
					// Not a C on ref
				}

			}
			
//			if(locus.getRecordAndPositions().size() > 0){
//				System.out.println(locus.getPosition() + " " + locus.getRecordAndPositions().size() + " " +
//						locus.getRecordAndPositions().get(0).getRecord().getFlags());
//			}
		}
	}
	
}
