package readWriteBAMUtils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import java.io.File;
import java.io.IOException;
import java.io.UnsupportedEncodingException;

/**
 * Convenience functions to read and write sam files supporting also to and from stdout.
 * TODO: Allow setting of SAMProgramRecords info.
 * @author berald01
 *
 */
public class ReadWriteBAMUtils {
	
	/**
	 * Open sam or bam file given filename or "-" to read from stdin.
	 * @param insam
	 * @param validationStringency See options in ValidationStringency. Typically one of:
	 * 		  ValidationStringency.SILENT, ValidationStringency.LENIENT, ValidationStringency.STRICT 
	 * @return
	 */
	public static SamReader reader(String insam, ValidationStringency validationStringency){
		SamReaderFactory sf = SamReaderFactory.
				makeDefault().
				validationStringency(validationStringency);
		SamReader sam= null;
		if(insam.equals("-")){
			SamInputResource resource= SamInputResource.of(System.in);
			sam= sf.open(resource);			
		} else {
			sam= sf.open(new File(insam));
		}
		return sam;
	}

	/**
	 * Write sam or bam file, depending on extension or to stdout.
	 * @param outsam Output filename, extension determines type. 
	 *     Use - or -.bam to write BAM to stdout or -.sam for SAM to stdout.  
	 * @param outHeader
	 * @return
	 * @throws UnsupportedEncodingException   
	 */
	public static SAMFileWriter writer(String outsam, SAMFileHeader outHeader) throws UnsupportedEncodingException {		
		SAMFileWriter outbam;
		outsam= outsam.trim();
		boolean presorted= false;
		if(outsam.equals("-") || outsam.equals("-.bam")){
			// TODO: Replace System.out as it is slow. 
			outbam= new SAMFileWriterFactory().makeBAMWriter(outHeader, presorted, System.out);
		} else if(outsam.equals("-.sam")){
			// TODO: Replace System.out as it is slow. 
			outbam= new SAMFileWriterFactory().makeSAMWriter(outHeader, presorted, System.out);
		} else {
			outbam= new SAMFileWriterFactory().
					makeSAMOrBAMWriter(outHeader, presorted, new File(outsam));
		}
		return outbam;
	}
	
	/**
	 * Open sam or bam file for writing taking header from template
	 * @param outsam File name to write to. Use - fro stdout bam or -.sam for stdout sam.
	 * @param template Filename of template sam or bam. Use header from this file.
	 * @return
	 * @throws IOException 
	 */
	public static SAMFileWriter writer(String outsam, String template) throws IOException{
		SamReaderFactory sf = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT);
		SamReader sam= sf.open(new File(template));
		SAMFileHeader header= sam.getFileHeader();		
		SAMFileWriter outbam= writer(outsam, header);
		sam.close();
		return outbam;
	}
}
