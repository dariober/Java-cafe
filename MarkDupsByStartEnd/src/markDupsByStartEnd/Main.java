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
		BufferedWriter br= new BufferedWriter(new FileWriter(tmp));
		
		// Read through sam file
		for(SAMRecord rec : sam){
			
			// Write out the reads unchanged that contain any one of these flags:
			if(rec.getReadPairedFlag() || 
			   rec.getReadUnmappedFlag() || 
			   rec.getReadFailsVendorQualityCheckFlag() ||
			   rec.getSupplementaryAlignmentFlag()){
				outbam.addAlignment(rec);
				nRecsSkipped++;
			} else {
				br.write(Utils.samRecordToTabLine(rec));
			}
		}
		br.close();
		// Sort tab separated file to bring together duplicates.
		InputStream sortedTabFile= Utils.sortTabAndGetOuput(tmp.getAbsolutePath()).getInputStream();
		
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
				dedupBlock= Arrays.copyOfRange(line, 0, 5);
				nRecsNonDups++;
			}
			outbam.addAlignment(rec);
		}
		
		outbam.close();
		System.err.println("N. records skipped\t" + nRecsSkipped);
		System.err.println("N. duplicates\t" + nRecsDups);
		System.err.println("N. non duplicates\t" + nRecsNonDups);
		System.exit(0);
	}

}
