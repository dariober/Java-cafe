package markDupsByStartEnd;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.zip.GZIPOutputStream;

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
import htsjdk.samtools.cram.encoding.writer.Writer;

public class Main {

	public static void main(String[] args) throws IOException, InterruptedException {
			
		int MAX_RECORDS_IN_RAM= 1000000;
		
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
					setMaxRecordsInRam(MAX_RECORDS_IN_RAM).
					makeSAMWriter(outhdr, false, System.out);	
		} else {
			outbam= new SAMFileWriterFactory().
					setMaxRecordsInRam(MAX_RECORDS_IN_RAM).
					makeSAMOrBAMWriter(outhdr, false, new File(outsam));
		}
		
		// Prepare tab delimited file that will be used to put together duplicate blocks
		// =============================================================================
		String tmpname= Utils.getTmpFilename(insam); 
		File tmp = new File(tmpname); // File.createTempFile("markDupsByStartEnd.", ".tmp.txt");
		tmp.deleteOnExit();
		BufferedWriter br= new BufferedWriter(new FileWriter(tmp));
		
		// System.err.println("Writing to\n" + tmp.getAbsolutePath()); // DEBUGGING
		
		// Read through sam file
		List<SAMRecordExt> lst= new ArrayList<SAMRecordExt>(); // STUB
		List<String> lstTmpFilenames= new ArrayList<String>(); //STUB
		int n= 0; //STUB
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
				br.write(Utils.samRecordToTabLine(rec, ignoreReadGroup));
				
				/* STUB: Accumulate records until you hit a certain limit (5M recs?)
				 * Then sort, write to tmp file. 
				 * See if ArrayList can be saved as ser and then read back one item at a time
				 * as you would for file.
				lst.add(new SAMRecordExt(rec, ignoreReadGroup));
				if(lst.size() > 5000000){
					Collections.sort(lst);
					OutputStreamWriter writer = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(tmpname + n)), "UTF-8");
					// BufferedWriter writer = new BufferedWriter(new FileWriter(tmpname + n));
					lstTmpFilenames.add(tmpname + n);
					for(SAMRecordExt x : lst){
						writer.write(x.toString());
					}
					writer.close();
					lst.clear();
					n++;
				} */
			}
			if(nRecsTot % 1000000 == 0){
				System.err.println("First pass: " + nRecsTot + " read.");
			}
		}
		br.close();
		// STUB: Code to write to file the last chunk of data in lst. Note that lst might be empty now. 
		// STUB: Code to read in parallel the files in lstTmpFilenames
		// STUB: Each top line from files converted to SAMREcordExt, put in a list and sorted
		// STUB: Process the stream of SAMREcordExt to pick mark duplicates,
		
		// Sort tab separated file to bring together duplicates.
		Process p= Utils.sortTabAndGetOuput(tmp.getAbsolutePath());
		InputStream sortedTabFile= p.getInputStream();

		// Mark duplicates
		// Thanks to the sorting above, the first read of each block is the
		// best one and is left unchanged. The following reads are marked.
		BufferedReader obr= new BufferedReader(new InputStreamReader(sortedTabFile));
		String str;
		String[] dedupBlock= null; 
		while ((str = obr.readLine()) != null) {
			String[] line= str.split("\t");
			String[] array= Arrays.copyOfRange(line, Utils.OFFSET, line.length);
			SAMRecord rec= Utils.arrayToSAMRecord(array, outhdr);
			String[] currentBlock= Arrays.copyOfRange(line, 0, 5);
			if(Arrays.deepEquals(currentBlock, dedupBlock)){
				rec.setDuplicateReadFlag(true);
				nRecsDups++;
			} else {
				rec.setDuplicateReadFlag(false); // <- NB: If read was marked dup. Now it is unmarked 
				dedupBlock= Arrays.copyOfRange(line, 0, 5);
				nRecsNonDups++;
			}
			outbam.addAlignment(rec);
			if((nRecsDups + nRecsNonDups) % 1000000 == 0){
				System.err.println(nRecsDups + nRecsNonDups + " records processed.");
			}
		}
		obr.close();
		outbam.close();
		if(p.exitValue() != 0){
			System.err.println("Sorting exited with error " + p.exitValue());
			BufferedReader ebr= new BufferedReader(new InputStreamReader(p.getErrorStream()));
			String x;
			while ((x = ebr.readLine()) != null) {
				System.err.println(x);
			}
			System.exit(1);
		}
		System.err.println("N. alignment total\t" + nRecsTot);
		System.err.println("N. alignment skipped\t" + nRecsSkipped);
		System.err.println("N. alignment duplicates\t" + nRecsDups);
		System.err.println("N. alignment non duplicates\t" + nRecsNonDups);
		System.exit(0);
	}

}
