package markDupsByStartEnd;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.zip.GZIPInputStream;
import net.sourceforge.argparse4j.inf.Namespace;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class Main {

	public static void main(String[] args) throws IOException, InterruptedException {
			
		int MAX_RECORDS_IN_RAM= 100000;
		
		/* Start parsing arguments */
		Namespace opts= ArgParse.argParse(args);
		ArgParse.validateArgs(opts);
	
		String insam= opts.getString("insam");
		String outsam= opts.getString("outsam");
		boolean unsortedOutput= opts.getBoolean("unsortedOutput");
		String vs= opts.getString("validationStringency");
		boolean ignoreReadGroup= opts.getBoolean("ignoreReadGroup");
		
		ValidationStringency validationStringency = null;
		if(vs.equals("SILENT")){
			validationStringency= ValidationStringency.SILENT;
		} else if(vs.equals("LENIENT")){
			validationStringency= ValidationStringency.LENIENT;
		} else if(vs.equals("STRICT")){
			validationStringency= ValidationStringency.STRICT;
		} else {
			System.err.println("Invalid string for validation stringency: " + vs);
			System.exit(1);
		}
			
		// Loggers
		// =======
		long nRecsSkipped= 0;
		long nRecsDups= 0;
		long nRecsNonDups= 0;
		long nRecsTot= 0;
			
		// Prepare to read input sam/bam
		SamReaderFactory sf = SamReaderFactory.makeDefault().validationStringency(validationStringency);
		SamReader sam= null;
		if(!insam.equals("-")){
			sam= sf.open(new File(insam));
		} else {
			SamInputResource resource= SamInputResource.of(System.in);
			sam= sf.open(resource);			
		}
		
		// Prepare output file
		// ===================
		SortOrder so= SortOrder.coordinate;
		if(unsortedOutput){
			so= SortOrder.unsorted;
		}
		SAMFileHeader outhdr= sam.getFileHeader();
		outhdr.setSortOrder(so);
		SAMFileWriter outbam;
		if(outsam.equals("-")){
			outbam= new SAMFileWriterFactory().
//					setMaxRecordsInRam(MAX_RECORDS_IN_RAM).
					makeSAMWriter(outhdr, false, System.out);	
		} else {
			outbam= new SAMFileWriterFactory().
					setMaxRecordsInRam(MAX_RECORDS_IN_RAM).
					makeSAMOrBAMWriter(outhdr, false, new File(outsam));
		}
		
		// Read through sam file
		// =====================
		List<SAMRecordExt> lst= new ArrayList<SAMRecordExt>();
		List<String> lstTmpFilenames= new ArrayList<String>();
		for(SAMRecord rec : sam){
									
			nRecsTot++;
			// Write out the reads unchanged that contain any one of these flags:
			if(rec.getReadPairedFlag() || 
			   rec.getReadUnmappedFlag() || 
			   rec.getReadFailsVendorQualityCheckFlag() ||
			   rec.getSupplementaryAlignmentFlag()){
				outbam.addAlignment(rec);
				nRecsSkipped++;
			} else {
				/* Accumulate records until you hit a certain limit (5M recs?)
				 * Then sort, write to tmp file. 
				 * See if ArrayList can be saved as ser and then read back one item at a time
				 * as you would for file. */
				
				try{
					lst.add(new SAMRecordExt(rec, ignoreReadGroup));
				} catch (NullPointerException e) {		
					System.err.println("Error reading read #" + nRecsTot 
							+ "\nAre the RG tags in reads are consistent with the RG dictionary in the sam header?\n"
							+ "To ignore the RG information use the --ignoreReadGroup/-rg option.\n"
							+ "\nStack Trace:");
					e.printStackTrace();
					System.err.printf("\nRead was\n%s\n", rec.getSAMString());
					System.exit(1);
				}
				
				if(lst.size() >= 500000){
					Collections.sort(lst);
					String tmpFileName= Utils.writeListToGzipFile(lst, insam);
					lstTmpFilenames.add(tmpFileName);
					lst.clear();
					
				}
			}
			//if(nRecsTot % 1000000 == 0){
			//	System.err.println("First pass: " + nRecsTot + " read.");
			//}

		}
		//Write to file the last chunk of data in lst. 
		Collections.sort(lst);		
		String tmpFileName= Utils.writeListToGzipFile(lst, insam);
		lstTmpFilenames.add(tmpFileName);
		
		System.err.println("Merge sorting and marking duplicates");
		
		// Code to read in parallel the files in lstTmpFilenames
		// Open a buffered reader for each file.
		List<BufferedReader> brLst= new ArrayList<BufferedReader>();
		for( String f : lstTmpFilenames ){
			InputStream fileStream = new FileInputStream(f);
			InputStream gzipStream = new GZIPInputStream(fileStream);
			Reader decoder = new InputStreamReader(gzipStream, "UTF-8");
			BufferedReader buffered = new BufferedReader(decoder);
			brLst.add(buffered);
		}
		
		// * Read top line from each file, put each item in list topLst.
        // * Return the min element from topLst xmin= Collections.min(topLst), 
		// * Get index of xmin i= topLst.indexOf(xmin)
		// * Pick another element from the file i, put the element where xmin was topLst.set(i, newXmin)
		// * Stop when all the files return null.
		List<SAMRecordExt> topLst= new ArrayList<SAMRecordExt>();
		String[] dedupBlk= null;
		for( BufferedReader ibr : brLst ){
			String line= ibr.readLine();
			if(line == null) {
				brLst.remove(brLst);
			} else {
				SAMRecordExt srec= new SAMRecordExt(line, outbam.getFileHeader());
				topLst.add(srec);
			}
		}
		long i= 0;
		while(topLst.size() > 0){
			
			SAMRecordExt xRecMin= Collections.min(topLst);
			int minidx= topLst.indexOf(xRecMin);
			// Refill topLst
			String line= brLst.get(minidx).readLine();
			if(line == null){
				brLst.get(minidx).close();
				topLst.remove(xRecMin);
				brLst.remove(brLst.get(minidx));
			} else {
				topLst.set(minidx, new SAMRecordExt(line, outbam.getFileHeader()));
			}
			
			// Mark duplicates
			String[] currentBlock= xRecMin.getBlockPosition();
			if(Arrays.deepEquals(currentBlock, dedupBlk)){
				xRecMin.getSamRecord().setDuplicateReadFlag(true);
				nRecsDups++;
			} else {
				xRecMin.getSamRecord().setDuplicateReadFlag(false); // <- NB: If read was marked dup. Now it is unmarked 
				dedupBlk= xRecMin.getBlockPosition();
				nRecsNonDups++;
			}
			outbam.addAlignment(xRecMin.getSamRecord());
			i++;
			if(i % 1000000 == 0){
				System.err.println(i + " to output");
			}
		}
		outbam.close();

		System.err.printf("N. alignment total\t%s\t%.2f\n", nRecsTot, 100.0 * nRecsTot/nRecsTot);
		System.err.printf("N. alignment skipped\t%s\t%.2f\n", nRecsSkipped, 100.0 * nRecsSkipped/nRecsTot);
		System.err.printf("N. alignment duplicates\t%s\t%.2f\n", nRecsDups, 100.0 * nRecsDups/nRecsTot);
		System.err.printf("N. alignment non duplicates\t%s\t%.2f\n", nRecsNonDups, 100.0 * nRecsNonDups/nRecsTot);
		System.exit(0);
	}
}
