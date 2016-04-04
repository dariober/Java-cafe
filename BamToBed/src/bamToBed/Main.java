package bamToBed;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.Iterator;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import net.sourceforge.argparse4j.inf.Namespace;

public class Main {

	public static void main(String[] args) throws IOException {

		Namespace opts= ArgParse.argParse(args);
		SamReader samReader= Utils.reader(opts.getString("inbam"), ValidationStringency.SILENT);
		
		Iterator<SAMRecord> sam= null;
		if( opts.getString("chrom").isEmpty() ){
			sam= samReader.iterator();
		} else {
			sam= samReader.query(opts.getString("chrom"), opts.getInt("from"), opts.getInt("to"), false);
		}
		
		BufferedWriter out = new BufferedWriter(new OutputStreamWriter(System.out));
		while(sam.hasNext()){
			SAMRecord rec= sam.next();
			if(rec.getReadUnmappedFlag()){ // Consistent with bedtools 
				continue;
			}
			if( ! Utils.filterSamRecord(rec, opts.getInt("requiredFlag"), opts.getInt("filterFlag"), opts.getInt("mapq"))){
				out.write(Utils.SAMRecordToBed(rec) + "\n");
			}
		}
		out.flush();
		samReader.close();
	}

}
