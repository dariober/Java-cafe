package readWriteBAMUtils;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
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
		new File("test/test_data/deleteme.bam").delete();
		
		// Write SAM
		out= ReadWriteBAMUtils.writer("test/test_data/deleteme.sam", sam.getFileHeader());
		out.addAlignment(rec);
		out.close();
		new File("test/test_data/deleteme.sam").delete();
		
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
		new File("test/test_data/deleteme-template.sam").delete();
	}
	
	@Test
	public void canAddProgramGroupToHeader() throws IOException{
		// Get a test file header:
		SamReaderFactory sf = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.STRICT);
		SamReader sam= sf.open(new File("test/test_data/header.sam"));
		SAMFileHeader header= sam.getFileHeader();
		SAMFileHeader newHeader= ReadWriteBAMUtils.addPGtoFileHeader(header, "progID", "myProg", "0.0.0", "myProg.jar -H");
		
		assertEquals("myProg", newHeader.getProgramRecord("progID").getProgramName());
		assertEquals("0.0.0", newHeader.getProgramRecord("progID").getProgramVersion());
		assertEquals("myProg.jar -H", newHeader.getProgramRecord("progID").getCommandLine());		
	}

}
