package vmatchToSam;

import static org.junit.Assert.assertEquals;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.List;

import exceptions.InvalidVmatchRecordException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import net.sourceforge.argparse4j.inf.Namespace;
import utils.Utils;

public class Main {

	public static void main(String[] args) throws IOException, InvalidVmatchRecordException, ClassNotFoundException, SQLException {

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
		
		VmatchSqlite vsql= new VmatchSqlite();
				
		vsql.getConn().setAutoCommit(false); // This is important: By default each insert is committed 
											 // as it is executed, which is slow. Let's commit in bulk at the end instead.
		
		System.err.println("Writing to temporary sqlite database: " + vsql.getSqlitedb().getAbsolutePath());
		int nrec= 0;
		PreparedStatement stmtInsert= vsql.getConn().prepareStatement(vsql.getSqlInsert());
		while( br.ready() ){
			// First accumulate records in sqlite
			VmatchRecord vmr= new VmatchRecord(br);
			if(!vmr.isEmpty()){
				vsql.add(vmr, stmtInsert);
				stmtInsert.execute();
				nrec += 1;
			}
			if(nrec % 100000 == 0){
				System.err.println(nrec + " written to db");
			}
		}
		stmtInsert.close();
		System.err.println(nrec + " written to db");
		vsql.getConn().commit();
		br.close();
		
		// Now retrieve records from db 
		Statement stmt= vsql.getConn().createStatement();
		System.err.println("Creating index");
		stmt.execute("CREATE INDEX qname_idx ON " + VmatchRecord.SQLITE_VMATCH_TABLE + " (queryNameOrIndex)");
		System.err.println("Getting unique query names");
		ResultSet rs= stmt.executeQuery(
	        "SELECT DISTINCT queryNameOrIndex FROM " + VmatchRecord.SQLITE_VMATCH_TABLE + 
	        " ORDER BY queryNameOrIndex");
		
		String sql= "SELECT * FROM " + VmatchRecord.SQLITE_VMATCH_TABLE +
			    " WHERE queryNameOrIndex = ?" +
			    " ORDER BY alignmentScore DESC";
		PreparedStatement qstmt= vsql.getConn().prepareStatement(sql);
		while(rs.next()){
			String qname= rs.getString("queryNameOrIndex");
			qstmt.setString(1, qname);
			ResultSet querySet= qstmt.executeQuery();
			List<SAMRecord> samRecStack= new ArrayList<SAMRecord>();
			while(querySet.next()){
				VmatchRecord vmr= new VmatchRecord(querySet);
				samRecStack.add(vmr.getSAMRecord());
			}
			Utils.setSecondaryAlignments(samRecStack);
			Utils.setMappingQualities(samRecStack);
			for(SAMRecord rec : samRecStack){
				System.out.print(rec.getSAMString());
			}
		}
		stmt.close();
	}

}
