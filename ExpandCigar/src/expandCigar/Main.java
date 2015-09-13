package expandCigar;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;

import java.io.IOException;
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
		
		/* Input */
		SamReader in= ReadWriteBAMUtils.reader(insam, ValidationStringency.SILENT);

		/* Set program group */
		SAMFileHeader header= in.getFileHeader();
		StringBuilder cmdLine= new StringBuilder();
		cmdLine.append(ArgParse.PROG_NAME);
		for(String x : args){
			cmdLine.append(" ");
			cmdLine.append(x);
		}
		header= ReadWriteBAMUtils.addPGtoFileHeader(header, 
				ArgParse.PROG_NAME, 
				ArgParse.PROG_NAME, 
				ArgParse.VERSION, 
				cmdLine.toString());				
		/* Output */
		SAMFileWriter out= ReadWriteBAMUtils.writer(outsam, header);
		
		/* Processing */
		for(SAMRecord rec : in){
			String xcigar= Utils.expandCigarMtoX(rec);
			if(oldCigarTag != null){
				rec.setAttribute(oldCigarTag, rec.getCigarString());
			}
			rec.setCigarString(xcigar);
			out.addAlignment(rec);
		}
		in.close();
		out.close();
	}

}
