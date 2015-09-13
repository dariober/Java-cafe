package expandCigar;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;

import java.io.IOException;
import java.util.Arrays;

import net.sourceforge.argparse4j.inf.Namespace;
import readWriteBAMUtils.ReadWriteBAMUtils;
import expandCigar.ArgParse;

public class Main {

	public static void main(String[] args) throws IOException {
		
		/* Parsing arguments */
		Namespace opts= ArgParse.argParse(args);
		ArgParse.validateArgs(opts);
	
		String insam= opts.getString("insam");
		String outsam= opts.getString("outsam");
		String oldCigarTag= opts.getString("oldCigarTag");
		
		/* Input/Output */
		SamReader in= ReadWriteBAMUtils.reader(insam, ValidationStringency.SILENT);
				
		SAMFileWriter out= ReadWriteBAMUtils.writer(outsam, in.getFileHeader());

		/* Set program group */
		
		// Prepare program group:
		SAMProgramRecord programRecord = new SAMProgramRecord(ArgParse.PROG_NAME);
		programRecord.setProgramName(ArgParse.PROG_NAME);
		programRecord.setCommandLine(Arrays.toString(args));
		programRecord.setProgramVersion(ArgParse.VERSION);
		// programRecord.setAttribute("XYZ", "Additional key-value pair");
		
		// Assign program group to header:
		SAMFileHeader header= out.getFileHeader();
		header.setProgramRecords(Arrays.asList(programRecord));
		
		/* Processing */
		for(SAMRecord rec : in){
			System.out.print(rec.getSAMString());
			String xcigar= Utils.expandCigarMtoX(rec);
			if(oldCigarTag != null){
				rec.setAttribute(oldCigarTag, rec.getCigarString());
			}
			rec.setCigarString(xcigar);
			System.out.print(rec.getSAMString());
			out.addAlignment(rec);
			
			break;
		}
		in.close();
		out.close();
	}

}
