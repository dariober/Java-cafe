package vmatchToSam;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;

import exceptions.InvalidVmatchRecordException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import net.sourceforge.argparse4j.inf.Namespace;
import utils.Utils;

public class Main {

	public static void main(String[] args) throws IOException, InvalidVmatchRecordException {

		Namespace opts= ArgParse.argParse(args);
		
		BufferedReader br;
		if(opts.getString("invmatch").equals("-") || opts.getString("invmatch").isEmpty()){
			br= new BufferedReader(new InputStreamReader(System.in));
		} else {
			br= new BufferedReader(new FileReader(new File(opts.getString("invmatch"))));
		}
		
		SAMFileHeader fh = null;
		if(!opts.getString("fai").isEmpty()){
			fh= Utils.makeSAMHeaderFromFastaIndex(opts.getString("fai"));
		} else if(!opts.getString("fasta").isEmpty()){
			fh= Utils.makeSAMHeaderFromFasta(opts.getString("fasta"));
		}
		if(fh != null){
			// All this is just to print the sam header to stdout.
			// since getTextHeader() returns null if called directly from object!!
			// First write the header to a tmp file, then read it back and print to stdout.
			SAMProgramRecord pr= new SAMProgramRecord(ArgParse.PROG_NAME);
			pr.setProgramVersion(ArgParse.VERSION);
			fh.addProgramRecord(pr);
			
			SAMReadGroupRecord rg= new SAMReadGroupRecord(VmatchRecord.DEFAULT_READ_GROUP);
			rg.setSample("NA");
			fh.addReadGroup(rg);
			
			File temp = File.createTempFile("vmatchToSam", ".tmp.sam");
			temp.deleteOnExit();
			SAMFileWriter tmpHdr= new SAMFileWriterFactory().makeSAMWriter(fh, true, temp);
			tmpHdr.close();
			SamReaderFactory sf = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT);
			SamReader sam= sf.open(temp);
			System.out.print(sam.getFileHeader().getTextHeader());
			sam.close();
		}
		
		while( br.ready() ){
			VmatchRecord vmr= new VmatchRecord(br);
			System.out.print(vmr.getSAMRecord().getSAMString());
		}
		br.close();
		
	}

}
