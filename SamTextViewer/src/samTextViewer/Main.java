package samTextViewer;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import net.sourceforge.argparse4j.inf.Namespace;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import readWriteBAMUtils.ReadWriteBAMUtils;

public class Main {

	private static String prettySeqPrinter(byte[] faSeq, boolean noFormat){
		
		String faSeqStr= "";
		if(!noFormat){
			for(int i= 0; i < faSeq.length; i++){
				char base= (char)faSeq[i];
				// For colour scheme see http://www.umass.edu/molvis/tutorials/dna/atgc.htm
				if(base == 'A' || base == 'a'){
					faSeqStr += "\033[34m" + base + "\033[0m";
				} else if(base == 'C' || base == 'c') {
					faSeqStr += "\033[31m" + base + "\033[0m";
				} else if(base == 'G' || base == 'g') {
					faSeqStr += "\033[32m" + base + "\033[0m";
				} else if(base == 'T' || base == 't') {
					faSeqStr += "\033[33m" + base + "\033[0m";
				} else {
					faSeqStr += base;
				} 
			}
		} else {
			faSeqStr= new String(faSeq);
		}
		return faSeqStr;
	}
	
	private static String ruler(int from, int to, int by){
		int current= from;
		String numberLine= "";
		while(current < to){
			if(current % by == 1){
				numberLine= numberLine + current;
				current += String.valueOf(current).length();
			} else {
				numberLine= numberLine + " ";
				current++;
			}
		}
		return numberLine;
	}

	/** Read all records in input bam intersecting the given coordinates, then 
	 * arrange them to fill up the window.
	 * @param bam
	 * @param chrom
	 * @param from
	 * @param to
	 * @param faSeq
	 * @return
	 */
	private static List<List<TextRead>> readAndStackSAMRecords(
			String bam, GenomicCoords gc, byte[] faSeq, int f_incl, int F_excl, int mapq, boolean bs) {

		SamReader samReader= ReadWriteBAMUtils.reader(bam, ValidationStringency.SILENT);
		Iterator<SAMRecord> sam= samReader.query(gc.getChrom(), gc.getFrom(), gc.getTo(), false);
		List<TextRead> textReads= new ArrayList<TextRead>();
		TextWindow textWindow= new TextWindow(gc.getFrom(), gc.getTo());
		while(sam.hasNext()){
			// First accumulate all reads in window
			SAMRecord rec= sam.next();
			if((rec.getFlags() & f_incl) == f_incl && (rec.getFlags() & F_excl) == 0 && rec.getMappingQuality() >= mapq){
				TextRead tr= new TextRead(rec, textWindow);
				if(bs){
					tr.setTextRead(tr.convertTextToBS(faSeq));
				} else {
					tr.setTextRead(tr.getTextChars(faSeq));
				}
				textReads.add(tr);				
			}
		}
		// Arrange reads nicely 
		List<List<TextRead>> stackReads= Utils.stackReads(textReads);
		return stackReads;
	}

	/** Print stack of reads
	 * @param stackReads Reads accumulated and arranged, typically by readAndStackSAMRecords();
	 * @param nReadsOut Max number of lines of output
	 * @param fmt Use formatted printing?
	 */
	private static void printer(List<List<TextRead>> stackReads, int maxLines, boolean noFormat) {
		int i= 0;
		while((maxLines < 0 || i < maxLines) && i < stackReads.size()){
			List<TextRead> line= stackReads.get(i);
			String outline= Utils.linePrinter(line, noFormat);
			System.out.println(outline);
			i++;
		}
	}
	/* ------------------- M A I N ------------------- */
	public static void main(String[] args) throws IOException {
		
		/* Start parsing arguments * 
		 * *** If you change something here change also in console input ***/
		Namespace opts= ArgParse.argParse(args);
		ArgParse.validateArgs(opts);
		
		List<String> insam= opts.getList("insam");
		String region= opts.getString("region");
		int windowSize= opts.getInt("windowSize") - 1; // Because end is included
		String fasta= opts.getString("fasta");
		int maxLines= opts.getInt("maxLines");
		int f_incl= opts.getInt("f");
		int F_excl= opts.getInt("F");
		int mapq= opts.getInt("mapq");
		boolean bs= opts.getBoolean("BSseq");
		boolean noFormat= opts.getBoolean("noFormat");
		
		if(fasta == null && bs == true){
			System.err.println("Warning:\n"
					+ "Bisulfite mode enabled without fasta reference. Methylated bases will not be shown.");
		}
				
		/* Prepare genomics coordinates to fetch*/
		GenomicCoords gc= null;
		if(region == null){ // Go to start of bam file
			gc= Utils.getStartCoordsOfBAM(insam.get(0));
			gc.setTo(gc.getFrom() + windowSize);
		} else {
			gc= GenomicCoords.goToRegion(region, insam.get(0), windowSize);
		}
		
		// -----------------------------------------
		while(true){
			byte[] faSeq= null;
			if(fasta != null){
				IndexedFastaSequenceFile faSeqFile = null;
				try {
					faSeqFile = new IndexedFastaSequenceFile(new File(fasta));
					try{
						faSeq= faSeqFile.getSubsequenceAt(gc.getChrom(), gc.getFrom(), gc.getTo()).getBases();
					} catch (NullPointerException e){
						System.err.println("Cannot fetch sequence " + gc.toString());
						e.printStackTrace();
					}
					faSeqFile.close();
				} catch (FileNotFoundException e) {
					e.printStackTrace();
				}
			} else {
				faSeq= new byte[gc.getTo() - gc.getFrom() + 1];
				for(int i= 0; i < faSeq.length; i++){
					faSeq[i]= 'N';
				}
			}
	
			String prettySeq= prettySeqPrinter(faSeq, noFormat);
			String prettyRuler= ruler(gc.getFrom(), gc.getTo(), 10);
			
			System.out.println(prettyRuler);
			System.out.println(prettySeq);
			for(String sam : insam){
				/*Always filter out read unmapped*/
				if((F_excl & 4) != 4){
					F_excl += 4;
				}
				List<List<TextRead>> stackReads= readAndStackSAMRecords(sam, gc, faSeq, f_incl, F_excl, mapq, bs);
				String header= gc.toString() + " - " + sam;
				if(!noFormat){
					header= "\033[0;34m" + header + "\033[0m";
				}
				System.out.println(header);
				printer(stackReads, maxLines, noFormat);
			}
			System.out.println(prettySeq);
			System.out.println(prettyRuler);
			
			/* Interactive input */
			System.err.print("\n[f]orward or [b]ack or move with -/+[int] e.g. -10k;\nUse short opts e.g. -F 16 -r <chr>:[from]; [q]uit: ");
			BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
			String rawInput= br.readLine().trim();
			/* Parse args */
			if(rawInput.equals("q")){
				break;
			}
			if(rawInput.equals("f")  
				|| rawInput.equals("b") 
				|| rawInput.matches("^\\-{0,1}\\d+.*") 
				|| rawInput.matches("^\\+{0,1}\\d+.*")){ // No cmd line args either f/b ops or ints
				String newRegion= Utils.parseConsoleInput(rawInput, gc).trim();
				gc= GenomicCoords.goToRegion(newRegion, insam.get(0), windowSize);				
			} else {
				List<String> clArgs= Arrays.asList(rawInput.split("\\s+"));
				if(clArgs.indexOf("-r") != -1){
					int i= clArgs.indexOf("-r") + 1;
					gc= GenomicCoords.goToRegion(clArgs.get(i), insam.get(0), windowSize);
				}
				if(clArgs.indexOf("-f") != -1){
					int i= clArgs.indexOf("-f") + 1;
					f_incl= Integer.parseInt(clArgs.get(i));
				}
				if(clArgs.indexOf("-F") != -1){
					int i= clArgs.indexOf("-F") + 1;
					F_excl= Integer.parseInt(clArgs.get(i));
				}				
				if(clArgs.indexOf("-q") != -1){
					int i= clArgs.indexOf("-q") + 1;
					mapq= Integer.parseInt(clArgs.get(i));
				}				
				if(clArgs.indexOf("-m") != -1){
					int i= clArgs.indexOf("-m") + 1;
					maxLines= Integer.parseInt(clArgs.get(i));
				}				
			}
			
			
		}
	}
}
