package readWriteBAMUtils;

import java.io.File;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.util.Arrays;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import org.junit.Test;

import readWriteBAMUtils.ReadWriteBAMUtils;

public class ReadWriteBAMUtilsTest {

	@Test
	public void canOpenSamOrBamFileByNameForReading() throws IOException {
		
		// TODO: Should test for reading stdin.
		SamReader sam= ReadWriteBAMUtils.reader("test/test_data/chr18.sam", ValidationStringency.STRICT);
		sam.iterator().next().getSAMString();
		sam.close();
		
		// Example usage to read through file.
		SamReader bam= ReadWriteBAMUtils.reader("test/test_data/chr18.bam", ValidationStringency.SILENT);
		for(SAMRecord rec : bam){
			rec.getSAMString();
		}
		bam.close();		
	}
	
	@Test
	public void canWriteBAMorSAMtoGivenFilename() throws UnsupportedEncodingException{
		SamReaderFactory sf = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.STRICT);
		SamReader sam= sf.open(new File("test/test_data/chr18.bam"));
		
		SAMRecord rec= sam.iterator().next();
		SAMFileWriter out;
		// Write BAM
		out= ReadWriteBAMUtils.writer("test/test_data/deleteme.bam", sam.getFileHeader());
		out.addAlignment(rec);
		out.close();
		
		// Write SAM
		out= ReadWriteBAMUtils.writer("test/test_data/deleteme.sam", sam.getFileHeader());
		out.addAlignment(rec);
		out.close();
		
		// Write SAM to stdout
		out= ReadWriteBAMUtils.writer("-.sam", sam.getFileHeader());
		out.addAlignment(rec);
		// out.close();

	}

	@Test
	public void canWriteBAMorSAMUsingTemplate() throws IOException{
		SAMFileWriter out= 
				ReadWriteBAMUtils.writer("test/test_data/deleteme-template.sam", "test/test_data/chr18.bam");
		out.close();
	}
	
	@Test
	public void exampleSettingProgramRecord() throws UnsupportedEncodingException{

		// Get a test file header:
		SamReaderFactory sf = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.STRICT);
		SamReader sam= sf.open(new File("test/test_data/header.sam"));
		SAMFileHeader header= sam.getFileHeader();
		
		// Prepare program group:
		SAMProgramRecord programRecord = new SAMProgramRecord("testerID");
		programRecord.setProgramName("tester");
		programRecord.setCommandLine("command line string here");
		programRecord.setProgramVersion("0.0.0");
		programRecord.setAttribute("XYZ", "Additional key-value pair");
		
		// Assign program group to header:
		header.setProgramRecords(Arrays.asList(programRecord));
		
		// Write to file and back
		SAMFileWriter out=
				ReadWriteBAMUtils.writer("test/test_data/tmp.sam", header);
		out.close();
		
		SamReader samhdr= ReadWriteBAMUtils.reader("test/test_data/tmp.sam", ValidationStringency.STRICT);
		System.out.println(samhdr.getFileHeader().getTextHeader());
	} 
}
