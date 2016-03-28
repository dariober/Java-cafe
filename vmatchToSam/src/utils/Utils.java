package utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;

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
}
