package markDupsByStartEnd;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Arrays;

import net.sourceforge.argparse4j.inf.Namespace;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

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
		
		// Prepare to read input sam/bam
		SamReader sam= SamReaderFactory.
				makeDefault().
				validationStringency(validationStringency).
				open(new File(insam));
		
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
		File tmp = File.createTempFile("markDupsByStartEnd.", ".tmp.txt");
		tmp.deleteOnExit();
		System.err.println("Writing to\n" + tmp.getAbsolutePath());
		BufferedWriter br= new BufferedWriter(new FileWriter(tmp));
		
		// Read through sam file
		int naln= 0;
		for(SAMRecord rec : sam){
			naln++;
			// Write out the reads unchanged that contain any one of these flags:
			if(rec.getReadPairedFlag() || 
			   rec.getReadUnmappedFlag() || 
			   rec.getReadFailsVendorQualityCheckFlag() ||
			   rec.getSupplementaryAlignmentFlag()){
				outbam.addAlignment(rec);
				nRecsSkipped++;
			} else {
				br.write(Utils.samRecordToTabLine(rec, ignoreReadGroup));
			}
			if(naln % 1000000 == 0){
				System.err.println(naln);
			}
		}
		br.close();
		System.err.println("N. records skipped\t" + nRecsSkipped);

		System.err.println("File size: " + new File(tmp.getAbsolutePath()).length());
		
		// Sort tab separated file to bring together duplicates.
		Process p= Utils.sortTabAndGetOuput(tmp.getAbsolutePath());
		InputStream sortedTabFile= p.getInputStream();

		// Mark duplicates
		// Thanks to the sorting above, the first read of each block is the
		// best one and is left unchanged. The following reads are marked.
		BufferedReader obr= new BufferedReader(new InputStreamReader(sortedTabFile));
		String str;
		String[] dedupBlock= null;
		int OFFSET= 8; // Index where sam record fields start. 0 based. 
		while ((str = obr.readLine()) != null) {
			String[] line= str.split("\t");
			String[] array= Arrays.copyOfRange(line, OFFSET, line.length);
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
			if((nRecsDups + nRecsNonDups + nRecsSkipped) % 1000000 == 0){
				System.err.println(nRecsDups + nRecsNonDups + nRecsSkipped + " records processed.");
			}
		}
		int pex= p.exitValue();
		if(pex != 0){
			System.err.println("Sorting exited with error " + pex);
			p.getErrorStream();
			// TODO: Code to return the exit message
		}
		
		outbam.close();
		System.err.println("N. duplicates\t" + nRecsDups);
		System.err.println("N. non duplicates\t" + nRecsNonDups);
		System.exit(0);
	}

}
