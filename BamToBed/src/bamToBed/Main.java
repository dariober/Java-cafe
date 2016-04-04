package bamToBed;

import java.io.IOException;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import net.sourceforge.argparse4j.inf.Namespace;

public class Main {

	public static void main(String[] args) throws IOException {

		Namespace opts= ArgParse.argParse(args);
		SamReader sam= Utils.reader(opts.getString("inbam"), ValidationStringency.SILENT);
		
		for(SAMRecord rec : sam){
			if( (rec.getFlags() & opts.getInt("requiredFlag")) != opts.getInt("requiredFlag") ){
				continue;
			}
			if( (rec.getFlags() & opts.getInt("filterFlag")) != 0 ){
				continue;
			} 
			if (rec.getMappingQuality() < opts.getInt("mapq")){
				continue;
			} 
			System.out.print(rec.getSAMString()); // TODO: Use BufferedWriter
		}
		sam.close();
	}

}
