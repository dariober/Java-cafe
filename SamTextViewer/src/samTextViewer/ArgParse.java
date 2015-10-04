package samTextViewer;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.Arguments;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;

public class ArgParse {
	
	public static String PROG_NAME= "SamTextViewer";
	public static String VERSION= "0.1.0";
	
	/* Parse command line args */
	public static Namespace argParse(String[] args){
		ArgumentParser parser= ArgumentParsers
				.newArgumentParser(PROG_NAME)
				.defaultHelp(true)
				.version("${prog} " + VERSION)
				.description("DESCRIPTION\n"
+ "Text based viewer for sam/bam files, hopefully more convenient then samtools tview\n"
+ " and with support for bisulfite converted reads.\n"
+ "\n"
+ "See also home page at https://github.com/dariober/Java-cafe/tree/master/SamTextViewer"
+ "");	
		parser.addArgument("--insam", "-i")
			.type(String.class)
			.required(true)
			.nargs("+")
			.help("Input bam files, must be sorted and indexed.");
		
		parser.addArgument("--region", "-r")
			.type(String.class)
			.required(false)
			.help("Region to visualize and to collect reads from.\n"
					+ "Format 1-based: chrom:start-end or chrom:start or chrom.\n"
					+ "If not given the view will be at the start of the first read found.");

		parser.addArgument("--windowSize", "-w")
			.type(Integer.class)
			.setDefault(160)
			.help("Window size to display. Ignored if --region is in format chrom:start-end");
		
		parser.addArgument("--fasta", "-fa")
			.type(String.class)
			.help("Optional reference fasta reference file.\n"
					+ "If given, must be indexed, e.g. with `samtools faidx ref.fa`");

		parser.addArgument("--f", "-f")
			.type(Integer.class)
			.setDefault(0)
			.help("Required sam flags");
		
		parser.addArgument("--F", "-F")
			.type(Integer.class)
			.setDefault(0)
			.help("Filtering sam flags");

		parser.addArgument("--mapq", "-q")
			.type(Integer.class)
			.setDefault(0)
			.help("Minumum mapping quality for a read to be printed (same as samtools view -q)");
		
		parser.addArgument("--maxLines", "-m")
			.type(Integer.class)
			.setDefault(-1)
			.help("Maximum number of lines to print for each read track. No limit If <0.");

		parser.addArgument("--maxDepthLines", "-d")
			.type(Integer.class)
			.setDefault(10)
			.help("Maximum number of lines to print for each coverage.");

		
		parser.addArgument("--BSseq", "-bs")
			.action(Arguments.storeTrue())
			.help("Bisulphite mode: Mark bases as methylated (M/m) or unmethylated (U/u).");

		
		parser.addArgument("--noFormat", "-nf")
			.action(Arguments.storeTrue())
			.help("Do not format output with non ascii chars (colour, bold, etc.)");

		parser.addArgument("--nonInteractive", "-ni")
			.action(Arguments.storeFalse())
			.help("Non interactive mode: Exit after having processed cmd line args.");

		
		parser.addArgument("--version", "-v").action(Arguments.version());
		
		Namespace opts= null;
		try{
			opts= parser.parseArgs(args);
		}
		catch(ArgumentParserException e) {
			parser.handleError(e);
			System.exit(1);
		}		
		return(opts);
	}
	
	public static void validateArgs(Namespace opts){
		// Optional checks to validate arguments.
	}
}
