package softClipBam;

import java.io.IOException;
import java.util.List;

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
		
		/* Input/Output */
		SamReader in= ReadWriteBAMUtils.reader(insam, ValidationStringency.SILENT);
		SAMFileWriter out= ReadWriteBAMUtils.writer(outsam, insam);

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
