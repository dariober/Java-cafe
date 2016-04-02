package utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.Statement;
import java.util.List;

import exceptions.InvalidVmatchRecordException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import vmatchToSam.VmatchRecord;

public class Utils {

	/** Create SAM header from reference fasta sequences. If file fasta has index, use index. 
	 * @throws IOException 
	 * */
	public static SAMFileHeader makeSAMHeaderFromFasta(String fastaName) throws IOException {
		
		if(new File(fastaName + ".fai").exists()){
			// Parse index
			SAMFileHeader fh= new SAMFileHeader();
			fh= makeSAMHeaderFromFastaIndex(fastaName + ".fai");
			return fh;
		} else {
			// Parse fasta to get sequence names and lengths
			BufferedReader br= new BufferedReader(new FileReader(new File(fastaName)));
			String line;
			String fname = null;
			int seqlen= 0;
			boolean start= true;
			SAMFileHeader fh= new SAMFileHeader();
			while( (line= br.readLine()) != null ){
				line= line.trim();
				if(line.startsWith(">")){
					if(!start){
						SAMSequenceRecord ssr= new SAMSequenceRecord(fname, seqlen);
						fh.addSequence(ssr);
						seqlen= 0;
					}
					start= false;
					fname= refSeqNameForSam(line.substring(1));
				} else {
					seqlen += line.length();
				}
			}
			fh.addSequence(new SAMSequenceRecord(fname, seqlen));
			br.close();
			return fh;
		}
	}

	/** Create SAM file header from fasta index. Any tab separated file with
	 * first two columns seqName and seqLength will do. 
	 * @throws IOException 
	 * */
	public static SAMFileHeader makeSAMHeaderFromFastaIndex(String fastaIndex) throws IOException {
		
		BufferedReader br= new BufferedReader(new FileReader(new File(fastaIndex)));
		String line; 
		SAMFileHeader fh= new SAMFileHeader();
		while((line= br.readLine())  != null){
			String[] xline= line.trim().split("\t");
			String name= refSeqNameForSam(xline[0]);
			int seqlen= Integer.parseInt(xline[1]);
			SAMSequenceRecord ssr= new SAMSequenceRecord(name, seqlen);
			fh.addSequence(ssr);
		}
		br.close();
		return fh;

	}

	/** Return input reference sequence name edited to comply with sam specs */
	public static String refSeqNameForSam(String refSeqName) {
		
		String clean= refSeqName.trim().replaceAll("\\s.*$", "");
		return clean;
	}

	/**
	 * Mark as supplementary the alignments with 2nd best AS.
	 * */
	public static void setSecondaryAlignments(List<SAMRecord> samRecStack) {
		if(samRecStack.size() > 1){
			for(int i= 1; i < samRecStack.size(); i++){
				samRecStack.get(i).setSupplementaryAlignmentFlag(true);
			}
		}		
	}

	/**
	 * Assign mapping quality to multiple alignments of the same query */
	public static void setMappingQualities(List<SAMRecord> samRecStack) {
		if(samRecStack.size() == 1){
			samRecStack.get(0).setMappingQuality(30);
		} else if(samRecStack.size() > 1){
			int firstAlnScore= (int) samRecStack.get(0).getAttribute("AS");
			int secondAlnScore= (int) samRecStack.get(1).getAttribute("AS");
			samRecStack.get(0).setMappingQuality(firstAlnScore - secondAlnScore);
			for(int i= 1; i < samRecStack.size(); i++){
				samRecStack.get(i).setMappingQuality(0);
			}
		}		
	}
	
}
