package softClipBam;

import java.io.IOException;
import java.util.List;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import readWriteBAMUtils.ReadWriteBAMUtils;
import softClipBam.ArgParse;
import net.sourceforge.argparse4j.inf.Namespace;

public class Main {

	public static void main(String[] args) throws IOException {

		/* Start parsing arguments */
		Namespace opts= ArgParse.argParse(args);
		ArgParse.validateArgs(opts);
	
		String insam= opts.getString("insam");
		String outsam= opts.getString("outsam");
		List<Object> clipRead1= opts.getList("clipRead1");
		List<Object> clipRead2= opts.getList("clipRead2");
		
		/* Input */
		SamReader in= ReadWriteBAMUtils.reader(insam, ValidationStringency.SILENT);
		/* Add program group */
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

		/* Clipping */
		long i= 0;
		for(SAMRecord rec : in){
			Clipper.clip(rec, 
					(Integer) clipRead1.get(0), 
					(Integer) clipRead1.get(1), 
					(Integer) clipRead2.get(0), 
					(Integer) clipRead2.get(1));
			out.addAlignment(rec);
			i++;
			if( i % 10000000 == 0 ){
				System.err.println("Reads processed: " + i);
			}
		}
		System.err.println("Total reads: " + i);
		out.close();
		System.exit(0);
	}
}
