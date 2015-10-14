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
import java.util.Random;

import org.apache.commons.lang3.StringUtils;

import coverageViewer.CoverageViewer;
import methylationViewer.MethylLoci;
import net.sourceforge.argparse4j.inf.Namespace;
import filter.FlagToFilter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.AggregateFilter;
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
					faSeqStr += "\033[107;34m" + base + "\033[0m";
				} else if(base == 'C' || base == 'c') {
					faSeqStr += "\033[107;31m" + base + "\033[0m";
				} else if(base == 'G' || base == 'g') {
					faSeqStr += "\033[107;32m" + base + "\033[0m";
				} else if(base == 'T' || base == 't') {
					faSeqStr += "\033[107;33m" + base + "\033[0m";
				} else {
					faSeqStr += base;
				} 
			}
		} else {
			faSeqStr= new String(faSeq);
		}
		return faSeqStr;
	}
	
	/**
	 * Count reads in interval using the given filters.
	 * NB: Make sure you use the same filter as in readAndStackSAMRecords();
	 * @param bam
	 * @param gc
	 * @param filters List of filters to apply
	 * @return
	 */
	private static long countReadsInWindow(String bam, GenomicCoords gc, List<SamRecordFilter> filters) {

		long cnt= 0;
		
		SamReader samReader= ReadWriteBAMUtils.reader(bam, ValidationStringency.SILENT);
		Iterator<SAMRecord> sam= samReader.query(gc.getChrom(), gc.getFrom(), gc.getTo(), false);
		AggregateFilter aggregateFilter= new AggregateFilter(filters);
		while(sam.hasNext()){
			SAMRecord rec= sam.next();
			if( !aggregateFilter.filterOut(rec) ){
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
			String bam, GenomicCoords gc, byte[] faSeq, List<SamRecordFilter> filters, boolean bs, int maxReadsStack) {

		SamReader samReader= ReadWriteBAMUtils.reader(bam, ValidationStringency.SILENT);
		
		long cnt= countReadsInWindow(bam, gc, filters);
		float probSample= (float)maxReadsStack / cnt;
		
		Iterator<SAMRecord> sam= samReader.query(gc.getChrom(), gc.getFrom(), gc.getTo(), false);
		TextWindow textWindow= new TextWindow(gc.getFrom(), gc.getTo());
		List<TextRead> textReads= new ArrayList<TextRead>();
		Random rand = new Random();
		
		AggregateFilter aggregateFilter= new AggregateFilter(filters); 
		while(sam.hasNext() && textReads.size() < maxReadsStack){

			SAMRecord rec= sam.next();
			if( !aggregateFilter.filterOut(rec) ){
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
	
	private static String getMemoryStat(){
		float mem= (float) ((float)Runtime.getRuntime().totalMemory() / 1000000d);
		String memStats= "Mem use: " +  Math.round(mem * 10)/10 + " MB";
		return memStats;
	}
	/* ------------------- M A I N ------------------- */
	public static void main(String[] args) throws IOException {
				
		/* Start parsing arguments * 
		 * *** If you change something here change also in console input ***/
		Namespace opts= ArgParse.argParse(args);
		
		List<String> insam= opts.getList("insam");
		String region= opts.getString("region");
		int windowSize= opts.getInt("windowSize");
		String fasta= opts.getString("fasta");
		int maxLines= opts.getInt("maxLines");
		int maxDepthLines= opts.getInt("maxDepthLines");
		int maxMethylLines= opts.getInt("maxMethylLines");
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
			gc.setTo(gc.getFrom() + windowSize - 1);
		} else {
			gc= GenomicCoords.goToRegion(region, insam.get(0), windowSize);
		}
		
		// -----------------------------------------
		while(true){ // Each loop processes the user's input.

			/* Ruler_TO_BE_DEPRECTED maps genomic coordinates to screen coordinates */
			// Ruler_TO_BE_DEPRECTED ruler= new Ruler_TO_BE_DEPRECTED(gc.getFrom(), gc.getTo(), windowSize);
			
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
			
			for(int i= 0; i < insam.size(); i++){ // Iterate through each input file
				String sam= insam.get(i);
				
				/* Coverage and methylation track                         */
				/* Coverage will be used for coverage and/or methylation  */
				/* ****************************************************** */
				CoverageViewer cw= null;
				if(maxDepthLines != 0 || (bs && maxMethylLines != 0)){  
					// Prepare if either coverage track or methylation track is required 
					cw= new CoverageViewer(sam, gc, windowSize, filters);
					cw.getMappingToScreen(gc.getMapping(windowSize));
				}
				/* Prepare methylation track */
				/* ========================= */
				String methylProfile= ""; // This is all you need to print out the methylation track.
				if(bs && maxMethylLines != 0){
					MethylLoci ml= null; 
					if(faSeq == null){
						faSeq= prepareRefSeq(fasta, gc);
					}
					ml= new MethylLoci(cw, faSeq, gc.getFrom());
					if(doCompress){
						ml.compressCovergeViewer(windowSize);	
					}
					/* Prepare printable methylation profile */
					/* ------------------------------------- */
					double depthPerMethylchar= (double)Math.round((float) ml.getMaxDepth() / maxMethylLines * 10d) / 10d;
					depthPerMethylchar= (depthPerMethylchar < 1) ? 1 : depthPerMethylchar;
					methylProfile += "Methylation depth; each '*.' = " + depthPerMethylchar + "x\n";
					methylProfile += StringUtils.join(ml.getMethylProfileStrings(maxMethylLines, noFormat), "\n") + "\n";
				}
				/* Prepare printable coverage track */
				/* ================================ */				
				String coverageTrack= ""; // This is all you need to print out coverage track.
				if(maxDepthLines != 0){
					//if(doCompress){
					//	cw.compressCovergeViewer(windowSize);	
					//}
					List<String> depthStrings= cw.getProfileStrings(maxDepthLines);
					double maxDepth= cw.getMaxDepth();
					double depthPerDot= (double)Math.round((double) maxDepth / maxDepthLines * 10d) / 10d;
					depthPerDot= (depthPerDot < 1) ? 1 : depthPerDot;
					coverageTrack= StringUtils.join(depthStrings, "\n");
	
					/* Coverage track Header */
					String fname= new File(sam).getName();
					String header= fname+ "; Max read depth: " + Math.round(maxDepth * 10d)/10d + "; Each . = " + depthPerDot + "x";
					if(!noFormat){
						header= "\033[0;34m" + header + "\033[0m";
					}
					coverageTrack= header + "\n"+ coverageTrack + "\n";
				}
				/*  ------------------------------------------------------------- */	

				/* Reads */
				/* ***** */
				String stackReadsStr= "";
				if(maxLines != 0 && !doCompress){
					List<List<TextRead>> stackReads= readAndStackSAMRecords(sam, gc, faSeq, filters, bs, maxReadsStack);
					stackReadsStr= stackReadsToString(stackReads, maxLines, noFormat);	
				}
				
				/* And print out... */
				/* **************** */
				System.out.print(coverageTrack);
				System.out.print(methylProfile);
				System.out.print(stackReadsStr);
								
			} // End loop through files 
			
			/* Footers and interactive prompt */
			/* ****************************** */
			
			/* Sequence */
			/* ======== */
			if(!doCompress){
				System.out.println(prettySeq);
			}
			// Print ruler 
			System.out.println(gc.printableRuler(windowSize, RULER_BY));

			/* Footer with window specs and filters*/ 
			/* ==================================== */
			String footer= gc.toString() + "; each char= " + Math.round(gc.getBpPerScreenColumn(windowSize) * 10d)/10d + " bp; " 
					+ "Filters: -q " + mapq  + " -f " + f_incl + " -F " + F_excl
					+ "; " + getMemoryStat();
			if(!noFormat){
				System.out.println("\033[0;34m" + footer + "\033[0m; ");
			} else {
				System.out.println(footer);
			}
			/* Interactive input */
			/* ================= */
			if(!nonInteractive){
				break;
			}
			String cmdInput= "";
			while(true){
				System.err.print("[f]orward, [b]ack; [zi]/[zo] zoom in/out; move by -/+[int][k|m] bases; [h] for help and options: ");
				BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
				cmdInput= br.readLine().trim();
				if(cmdInput.equals("")){
					// Nothing to do re-ask what to do.
				} else if(cmdInput.equals("h")){
					String inline= "\nNavigation options\n\n"
							+ "[f]/[b]:   Small step forward/backward (1/10 of a window)\n"
							+ "[ff]/[bb]: Large step forward/backward (1 window)\n"
							+ "[zi]/[zo]: Zoom in / zoom out\n"
							+ "[:pos]:    Go to position <pos> on current chromosome\n" 
							+ "+/- <int>: Move forward/backward this <int> number of bases. Suffixes k and m are allowed. E.g. '-2m' moves 2 Mbp back\n"
							+ "[q]:       Quit\n"
							+ "[h]:       Show this help";
					// [f]orward, [b]ack; Zoom in/out with [zi]/[zo]; Jump to pos -/+[int][k|m]
					System.out.println(inline);
					System.out.println("\nFiltering and visualization options\n");
					System.out.println(ArgParse.getDocstrings());
				} else {
					break;
				}
			}
			/* Parse args */
			/* ---------- */
			if(cmdInput.equals("q")){
				break;
			}
			if(cmdInput.equals("f") 
				|| cmdInput.equals("b")
				|| cmdInput.equals("ff") 
				|| cmdInput.equals("bb")
				|| cmdInput.startsWith(":")
				|| cmdInput.matches("^\\-{0,1}\\d+.*") 
				|| cmdInput.matches("^\\+{0,1}\\d+.*")){ // No cmd line args either f/b ops or ints
				cmdInput= cmdInput.matches("^\\+.*") ? cmdInput.substring(1) : cmdInput;
				String newRegion= Utils.parseConsoleInput(cmdInput, gc).trim();
				gc= GenomicCoords.goToRegion(newRegion, insam.get(0), windowSize);
			} else if(cmdInput.equals("zo")){
				gc.zoomOut();
			} else if(cmdInput.equals("zi")){
				gc.zoomIn();
			} else {
				List<String> clArgs= Arrays.asList(cmdInput.split("\\s+"));
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
				if(clArgs.indexOf("-ml") != -1){
					int i= clArgs.indexOf("-ml") + 1;
					maxMethylLines= Integer.parseInt(clArgs.get(i));
				}
			} 
			System.out.println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
		} // End while loop keep going until quit or if no interactive input set
	}
}
