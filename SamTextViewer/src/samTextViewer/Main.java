package samTextViewer;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Random;

import org.apache.commons.lang3.StringUtils;

import coverageViewer.CoverageViewer;
import net.sourceforge.argparse4j.inf.Namespace;
import filter.FlagToFilter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.MappingQualityFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import readWriteBAMUtils.ReadWriteBAMUtils;

/**
 * @author berald01
 *
 */
public class Main {

	/** Print the genomic position every so many text chars */
	final static int RULER_BY= 10; 
//	/** Stack up to so many reads to print in read window */
//	final static int MAX_READS_IN_STACK= 2000;
	
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

	/**
	 * Count reads in interval using the given filters.
	 * NB: Make sure you use the same filter as in readAndStackSAMRecords();
	 * @param bam
	 * @param gc
	 * @param f_incl
	 * @param F_excl
	 * @param mapq
	 * @return
	 */
	private static long countReadsInWindow(String bam, GenomicCoords gc, int f_incl, int F_excl, int mapq) {

		long cnt= 0;
		
		SamReader samReader= ReadWriteBAMUtils.reader(bam, ValidationStringency.SILENT);
		Iterator<SAMRecord> sam= samReader.query(gc.getChrom(), gc.getFrom(), gc.getTo(), false);
		while(sam.hasNext()){
			SAMRecord rec= sam.next();
			if((rec.getFlags() & f_incl) == f_incl && (rec.getFlags() & F_excl) == 0 && rec.getMappingQuality() >= mapq){
				cnt++;
			}
		}
		try {
			samReader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return cnt;
	}

	/** Read all records in input bam intersecting the given coordinates, then 
	 * arrange them to fill up the window. If you change filters here change also 
	 * in the counter!
	 * @param bam
	 * @param chrom
	 * @param from
	 * @param to
	 * @param faSeq
	 * @return
	 */
	private static List<List<TextRead>> readAndStackSAMRecords(
			String bam, GenomicCoords gc, byte[] faSeq, int f_incl, int F_excl, int mapq, boolean bs, int maxReadsStack) {

		SamReader samReader= ReadWriteBAMUtils.reader(bam, ValidationStringency.SILENT);
		
		long cnt= countReadsInWindow(bam, gc, f_incl, F_excl, mapq);
		float probSample= (float)maxReadsStack / cnt;
		
		Iterator<SAMRecord> sam= samReader.query(gc.getChrom(), gc.getFrom(), gc.getTo(), false);
		TextWindow textWindow= new TextWindow(gc.getFrom(), gc.getTo());
		List<TextRead> textReads= new ArrayList<TextRead>();
		Random rand = new Random();
		while(sam.hasNext() && textReads.size() < maxReadsStack){

			SAMRecord rec= sam.next();
			if((rec.getFlags() & f_incl) == f_incl && (rec.getFlags() & F_excl) == 0 && rec.getMappingQuality() >= mapq){
				
				if(rand.nextFloat() < probSample){ // Downsampler
					TextRead tr= new TextRead(rec, textWindow);
					if(bs){
						tr.setTextRead(tr.convertTextToBS(faSeq));
					} else {
						tr.setTextRead(tr.getTextChars(faSeq));
					}
					textReads.add(tr);
				}
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
	private static String stackReadsToString(List<List<TextRead>> stackReads, int maxLines, boolean noFormat) {
		int i= 0;
		StringBuilder sb= new StringBuilder();
		while((maxLines < 0 || i < maxLines) && i < stackReads.size()){
			List<TextRead> line= stackReads.get(i);
			sb.append(Utils.linePrinter(line, noFormat)).append("\n");
			i++;
		}
		return sb.toString();
	}
	
	/**
	 * Get sequence as byte[] for the given genomic coords. Array of N if fasta is null.
	 * @param fasta
	 * @param gc
	 * @return
	 * @throws IOException
	 */
	private static byte[] prepareRefSeq(String fasta, GenomicCoords gc) throws IOException{

		byte[] faSeq= null;
		if(fasta != null){
			IndexedFastaSequenceFile faSeqFile = null;
			try {
				faSeqFile = new IndexedFastaSequenceFile(new File(fasta));
				try{
					System.out.println(gc);
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
		return faSeq;
	} 
	
	/* ------------------- M A I N ------------------- */
	public static void main(String[] args) throws IOException {
		
		/* Start parsing arguments * 
		 * *** If you change something here change also in console input ***/
		Namespace opts= ArgParse.argParse(args);
		ArgParse.validateArgs(opts);
		
		List<String> insam= opts.getList("insam");
		String region= opts.getString("region");
		final int windowSize= opts.getInt("windowSize");
		String fasta= opts.getString("fasta");
		int maxLines= opts.getInt("maxLines");
		int maxDepthLines= opts.getInt("maxDepthLines");
		final int maxReadsStack= opts.getInt("maxReadsStack");
		int f_incl= opts.getInt("f");
		int F_excl= opts.getInt("F");
		int mapq= opts.getInt("mapq");
		boolean bs= opts.getBoolean("BSseq");
		boolean noFormat= opts.getBoolean("noFormat");
		boolean nonInteractive= opts.getBoolean("nonInteractive");
		
		if((F_excl & 4) != 4){ // Always filter out read unmapped
			F_excl += 4;
		}
		
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

			/* Prepare filters */
			
			List<SamRecordFilter> filters= FlagToFilter.flagToFilterList(f_incl, F_excl); // new ArrayList<SamRecordFilter>();
			filters.add(new MappingQualityFilter(mapq));
			
			/* Need to compress? */
			boolean doCompress= false;
			if((gc.getTo() - gc.getFrom()) > windowSize){
				doCompress= true;
			} 
			
			/* Prepare reference seq if necessary */
			byte[] faSeq= null;
			String prettySeq= null;
			if(!doCompress){ // Print sequence only if compression is not needed. 
				faSeq= prepareRefSeq(fasta, gc);	
				prettySeq= prettySeqPrinter(faSeq, noFormat);
				System.out.println(prettySeq);
			}
			
			String prettyRuler= "";
			for(int i= 0; i < insam.size(); i++){
				String sam= insam.get(i);
				
				/* coverage track */
				int maxDepth= -1;
				double depthPerLine= -1;
				String depthTrack= "";
				if(maxDepthLines != 0){
					CoverageViewer cw= new CoverageViewer(sam, gc.getChrom(), gc.getFrom(), gc.getTo(), 
							windowSize, filters);
					if(doCompress){
						cw.compressCovergeViewer(windowSize);
						prettyRuler= cw.ruler(RULER_BY);
					}
					List<String> depthStrings= cw.getProfileStrings(maxDepthLines);
					maxDepth= cw.getMaxDepth();
					depthPerLine= (double)Math.round((float) maxDepth / maxDepthLines * 10d) / 10d;
					depthPerLine= (depthPerLine < 1) ? 1 : depthPerLine;
					depthTrack= StringUtils.join(depthStrings, "\n") + "\n";
				}

				/* Prepare and print ruler */
				if(!doCompress){
					prettyRuler= ruler(gc.getFrom(), gc.getTo(), RULER_BY);			
				}
				System.out.println(prettyRuler);
				
				/* Header */
				String fname= new File(sam).getName();
				String header= fname+ "; Max read depth: " + maxDepth + "; Each '*': " + depthPerLine + "x";
				if(!noFormat){
					header= "\033[0;34m" + header + "\033[0m";
				}
				
				/* Reads */
				String stackReadsStr= "";
				if(maxLines != 0 && !doCompress){
					List<List<TextRead>> stackReads= readAndStackSAMRecords(sam, gc, faSeq, f_incl, F_excl, mapq, bs, maxReadsStack);
					stackReadsStr= stackReadsToString(stackReads, maxLines, noFormat);	
				}
				/* And print out... */
				System.out.println(header);
				System.out.print(depthTrack);
				System.out.print(stackReadsStr);
			}
			if(!doCompress){
				System.out.println(prettySeq);
			}
			System.out.println(prettyRuler);
			
			/* Footer */ 
			String footer= gc.toString() + "; -q " + mapq  + " -f " + f_incl + " -F " + F_excl;
			if(!noFormat){
				System.out.println("\033[0;34m" + footer + "\033[0m; ");
			} else {
				System.out.println(footer);
			}
			/* Interactive input */
			if(!nonInteractive){
				break;
			}
			System.err.print("[f]orward, [b]ack; Zoom in/out with [zi]/[zo]; Jump to pos -/+[int][k|m] e.g. +1m or -10k;\nSet cmd line short opts e.g. -F 16 -r <chr>:[from]; [q]uit: ");
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
				rawInput= rawInput.matches("^\\+.*") ? rawInput.substring(1) : rawInput;
				String newRegion= Utils.parseConsoleInput(rawInput, gc).trim();
				gc= GenomicCoords.goToRegion(newRegion, insam.get(0), windowSize);
			} else if(rawInput.equals("zo")){
				gc.zoomOut();
			} else if(rawInput.equals("zi")){
				gc.zoomIn();
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
					if((F_excl & 4) != 4){ // Always filter out read unmapped
						F_excl += 4;
					}
				}				
				if(clArgs.indexOf("-q") != -1){
					int i= clArgs.indexOf("-q") + 1;
					mapq= Integer.parseInt(clArgs.get(i));
				}				
				if(clArgs.indexOf("-m") != -1){
					int i= clArgs.indexOf("-m") + 1;
					maxLines= Integer.parseInt(clArgs.get(i));
				}				
				if(clArgs.indexOf("-d") != -1){
					int i= clArgs.indexOf("-d") + 1;
					maxDepthLines= Integer.parseInt(clArgs.get(i));
				}				

			}	
		}
	}
}
