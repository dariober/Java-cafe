package expandCigar;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.Arguments;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;

public class ArgParse {
	
	public static String VERSION= "0.1.0";
	public static String PROG_NAME= "ExpandCigar";
	
	/* Parse command line args */
	public static Namespace argParse(String[] args){
		ArgumentParser parser= ArgumentParsers
				.newArgumentParser(PROG_NAME)
				.defaultHelp(true)
				.version("${prog} " + VERSION)
				.description("DESCRIPTION\n"
+ "Expand cigar string in inout bam file to convert M operator to X and =.\n"
+ "");	
		parser.addArgument("--insam", "-i")
			.type(String.class)
			.required(true)
			.help("Input sam or bam file. Use - to read from stdin.");
		
		parser.addArgument("--outsam", "-o")
			.type(String.class)
			.required(false)
			.setDefault("-")
			.help("Output file. Format will be sam or bam depending on extension.\n"
					+ "Use - to print BAM to stdout or '.sam' for sam to stdout.");
		
		parser.addArgument("--oldCigarTag", "-t")
			.type(String.class)
			.required(false)
			.help("Store the original cigar in this tag, e.g. 'XX'. If this tag exists it will be overwritten.");
		
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
