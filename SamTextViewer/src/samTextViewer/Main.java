package samTextViewer;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.text.StrTokenizer;

import exceptions.InvalidGenomicCoordsException;
import net.sourceforge.argparse4j.inf.Namespace;
import filter.FlagToFilter;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.filter.MappingQualityFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import jline.console.ConsoleReader;
import jline.console.completer.StringsCompleter;
import tracks.IntervalFeatureSet;
import tracks.TrackCoverage;
import tracks.TrackIntervalFeature;
import tracks.TrackMethylation;
import tracks.TrackReads;
import tracks.TrackWiggles;

/**
 * @author berald01
 *
 */
public class Main {
	
	private static String getMemoryStat(){
		float mem= (float) ((float)Runtime.getRuntime().totalMemory() / 1000000d);
		String memStats= "Mem: " +  Math.round(mem * 10)/10 + " MB";
		return memStats;
	}
	/* ------------------- M A I N ------------------- */
	public static void main(String[] args) throws IOException, InvalidGenomicCoordsException {
		
		/* Start parsing arguments * 
		 * *** If you change something here change also in console input ***/
		Namespace opts= ArgParse.argParse(args);
		
		List<String> insam= opts.getList("insam");
		String region= opts.getString("region");
		int windowSize= opts.getInt("windowSize");
		String fasta= opts.getString("fasta");
		boolean rpm= opts.getBoolean("rpm");
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
		boolean withReadName= false; // FIXME: Add to parser?
		
		if((F_excl & 4) != 4){ // Always filter out read unmapped
			F_excl += 4;
		}
			
		if(fasta == null && bs == true){
			System.err.println("Warning: Fasta reference not provided. Bisulfite mode will be disabled");
			bs= false;
		}
		
		/* Test input files exist */
		List<String> dropMe= new ArrayList<String>();
		for(String x : insam){
			if(!new File(x).exists()){
				System.err.println("\nWarning: Dropping file " + x + " as it does not exist.\n");
				dropMe.add(x);
			}
		}
		for(String x : dropMe){
			insam.remove(x);
		}
		
		if(insam.size() == 0 && fasta == null){
			System.err.println("\nNo files in input: Nothing to be done!\n");
			System.exit(1);
		}
		
		/* Prepare genomics coordinates to fetch */
		if(region.isEmpty()){
			for(String x : insam){
				try {
					region= Utils.initRegionFromFile(x);
					break;
				} catch(Exception e){
					System.err.println("Could not initilize from file " + x);
				}
			}
		}
		if(region.isEmpty() && fasta != null){ // Try to initilize from fasta
			IndexedFastaSequenceFile faSeqFile = new IndexedFastaSequenceFile(new File(fasta));
			region= faSeqFile.nextSequence().getName();
			faSeqFile.close();
		}
		GenomicCoordsHistory gch= new GenomicCoordsHistory();
		SAMSequenceDictionary samSeqDict = GenomicCoords.getSamSeqDictFromAnyFile(insam, fasta);
		gch.add(new GenomicCoords(region, samSeqDict, windowSize, fasta));
		
		/* Files to parse as features (bed, gtf, etc) */
		LinkedHashMap<String, IntervalFeatureSet> intervalFiles= Utils.createIntervalFeatureSets(insam); 
		
		ConsoleReader console = new ConsoleReader();
		for(String x : insam){
			console.addCompleter(new StringsCompleter(new File(x).getName()));
		}
		
		boolean printIntervalFeatures= false;
		while(true){ // Each loop processes the user's input.
			
			/* Prepare filters */
			List<SamRecordFilter> filters= FlagToFilter.flagToFilterList(f_incl, F_excl); // new ArrayList<SamRecordFilter>();
			filters.add(new MappingQualityFilter(mapq));
			
			for(int i= 0; i < insam.size(); i++){ // Iterate through each input file
				
				String sam= insam.get(i);
				String title= new File(sam).getName() + "; ";
				if(Utils.getFileTypeFromName(sam).equals("bam")){
					/* Coverage and methylation track */
					String printableTrackCoverage= "";
					String printableTrackMethylation= "";
					if(maxDepthLines > 0 || bs){
						TrackCoverage trackCoverage= new TrackCoverage(sam, gch.current(), filters, bs);
						if(maxDepthLines > 0){
							String printable= trackCoverage.printToScreen(gch.current().getMapping(), maxDepthLines, rpm);
							title += "Each . = " + Math.rint(trackCoverage.getScorePerDot()*100)/100 +  "; max depth: " + Math.rint(trackCoverage.getMaxDepth()*100)/100 + "x; ";
							printableTrackCoverage += printable + "\n";
						}
						if(bs && maxMethylLines > 0){
							TrackMethylation trackMethylation= new TrackMethylation(trackCoverage.getScreenLocusInfoList());
							printableTrackMethylation += trackMethylation.printToScreen(maxMethylLines, noFormat) + "\n";
							title += "BS[each */. = " + Math.rint(trackMethylation.getScorePerDot()*100)/100 + "x; " + 
									"max depth: " + Math.rint(trackMethylation.getMaxDepth()*100)/100 + "x" + "]; ";
						}
					}
					
					/* Reads */
					String printableTrackReads= "";
					if(maxLines != 0){
						TrackReads trackReads= new TrackReads(sam, gch.current(), filters, maxReadsStack);
						printableTrackReads += trackReads.printToScreen(maxLines, bs, noFormat, withReadName) + "\n";
					}								
					
					/* Print to screen */
					if(noFormat){
						System.out.println(title);
					} else {
						System.out.println("\033[0;34m" + title + "\033[0m");
					}
					System.out.print(printableTrackCoverage);
					System.out.print(printableTrackMethylation);
					System.out.print(printableTrackReads);
				} // End processing bam file
				if(intervalFiles.containsKey(sam)){
					TrackIntervalFeature tif= new TrackIntervalFeature(intervalFiles.get(sam), gch.current());
					String printableIntervalTrack= tif.printToScreen(noFormat);
					System.out.println(new File(sam).getName());
					System.out.println(printableIntervalTrack);
					if(printIntervalFeatures){
						System.out.println(tif.toString());
					}
				} 
				if(Utils.getFileTypeFromName(sam).equals("bigWig") 
						|| Utils.getFileTypeFromName(sam).equals("tdf")
						|| Utils.getFileTypeFromName(sam).equals("bedGraph")){
					TrackWiggles tw= new TrackWiggles(sam, gch.current());
					String printable= tw.printToScreen(maxDepthLines);
					//	title += "Each . = " + Math.rint(trackCoverage.getScorePerDot()*100)/100 +  "; max depth: " + Math.rint(trackCoverage.getMaxDepth()*100)/100 + "x; ";
					title += "Each . = " + Math.rint(tw.getScorePerDot()*1000)/1000 +  "; max score: " + Math.rint(tw.getMaxDepth()*1000)/1000 + "; ";
					if(noFormat){
						System.out.println(title);
					} else {
						System.out.println("\033[0;34m" + title + "\033[0m");
					}
					System.out.println(printable);
				}
			} // End loop through files 
			
			/* Footers and interactive prompt */
			/* ****************************** */
			System.out.print(gch.current().printableRefSeq(noFormat));
			System.out.println(gch.current().printableRuler(10));

			String footer= gch.current().toString() + "; " + Math.rint(gch.current().getBpPerScreenColumn() * 10d)/10d + " bp/char; " 
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
			while(cmdInput.isEmpty()){ // Keep asking for input until you get something valid
				console.setPrompt("[h] for help: ");
				cmdInput = console.readLine().trim();
				
				if(cmdInput.equals("h")){
					String inline= "\nNavigation options\n\n"
							+ "[f]/[b]\n        Small, step forward/backward 1/10 window\n"
							+ "[ff]/[bb]\n        Large step forward/backward 1/2 window\n"
							+ "[zi]/[zo]\n        Zoom in / zoom out\n"
							+ "[p]/[n]\n       Go to previous/next visited position\n"
							+ "[:pos]\n        Go to position <pos> on current chromosome\n" 
							+ "+/-<int>\n        Move forward/backward by <int> bases. Suffixes k and m allowed. E.g. -2m\n"
							+ "<filename> [:n]\n        Move to the next feature in <filename> on current chromosome\n"
							+ "<filename> <string> [:find]\n        Find the next record in file containing string. Use single quotes for strings containing spaces.\n"
							+ "[print]\n        Turn on/off the printing of bed/gtf features in current interval\n"
							+ "[rNameOn]/[rNameOff]\n        Show/Hide read names\n"
							+ "[history]\n        Show visited positions\n"
							+ "[q]\n        Quit\n"
							+ "[h]\n        Show this help";
					System.out.println(inline);
					System.out.println("\nCommand line options\n");
					System.out.println(ArgParse.getDocstrings());
					cmdInput= "";
					continue;
				} 
				/* Parse args */
				/* ---------- */
				try{
					if(cmdInput.equals("q")){
						System.exit(0);
					}
					
					if(cmdInput.equals("f") 
						|| cmdInput.equals("b")
						|| cmdInput.equals("ff") 
						|| cmdInput.equals("bb")
						|| cmdInput.matches("^:\\d+")
						|| cmdInput.matches("^\\-{0,1}\\d+.*") 
						|| cmdInput.matches("^\\+{0,1}\\d+.*")){ // No cmd line args either f/b ops or ints
						cmdInput= cmdInput.matches("^\\+.*") ? cmdInput.substring(1) : cmdInput;
						String newRegion= Utils.parseConsoleInput(cmdInput, gch.current()).trim();
						gch.add(new GenomicCoords(newRegion, samSeqDict, windowSize, fasta));
					} else if (cmdInput.equals("p")) {
						gch.previous();
					} else if (cmdInput.equals("n")) {
						gch.next();
					} else if(cmdInput.equals("zo")){
						GenomicCoords gc= (GenomicCoords)gch.current().clone();
						gc.zoomOut();
						gch.add(gc);
					} else if(cmdInput.equals("zi")){
						GenomicCoords gc= (GenomicCoords)gch.current().clone();
						gc.zoomIn();
						gch.add(gc);
					} else if(cmdInput.equals("history")){
						for(GenomicCoords x : gch.getHistory()){
							System.out.println(x);
						}
						cmdInput= "";
					} else if(cmdInput.toLowerCase().equals("print")){
						printIntervalFeatures= (printIntervalFeatures) ? false : true;
						System.err.println("Print interval features: " + printIntervalFeatures);
					} else if(cmdInput.toLowerCase().equals("rnameon")){
						withReadName= true;
					} else if(cmdInput.toLowerCase().equals("rnameoff")) {
						withReadName= false;
					} else if(cmdInput.endsWith(":n")){
						GenomicCoords gc= (GenomicCoords)gch.current().clone();
						gch.add(Utils.goToNextFeatureOnFile(cmdInput.replaceAll(":n$", "").trim(), gc, intervalFiles));
					} else if(cmdInput.endsWith(" :find")) {  
					
						StrTokenizer str= new StrTokenizer(cmdInput);
						str.setQuoteChar('\'');
						List<String> tokens= str.getTokenList();
						if(tokens.size() != 3){
							System.err.println("Error in :find subcommand. Expected 3 args got: " + cmdInput);
							cmdInput= "";
							continue;
						}
						GenomicCoords gc= (GenomicCoords)gch.current().clone();
						gch.add(Utils.findNextStringOnFile(tokens.get(1), tokens.get(0), gc, intervalFiles));
					
					} else { // Command line options from Argparse
						List<String> clArgs= Arrays.asList(cmdInput.split("\\s+"));
						if(clArgs.indexOf("-r") != -1){
							int i= clArgs.indexOf("-r") + 1;
							gch.add(new GenomicCoords(clArgs.get(i), samSeqDict, windowSize, fasta));
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
						if(clArgs.indexOf("-rpm") != -1){
							rpm= (rpm) ? false : true; // Invert rpm set.
						}
					} // END Command line options from Argparse
				} catch(Exception e){
					System.err.println("\nError processing input: " + cmdInput);
					System.err.println("\nStack trace:\n");
					e.printStackTrace();
					System.err.println("");
					cmdInput= "";
				}
			} // END while loop to get cmLine args
				System.out.println(StringUtils.repeat("~", gch.current().getUserWindowSize()));
		} // End while loop keep going until quit or if no interactive input set
	}
}
