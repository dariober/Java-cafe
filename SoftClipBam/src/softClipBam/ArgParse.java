package softClipBam;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.Arguments;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;

public class ArgParse {
	
	public static String VERSION= "0.1.0";
	
	/* Parse command line args */
	public static Namespace argParse(String[] args){
		ArgumentParser parser= ArgumentParsers
				.newArgumentParser("SoftClipBam")
				.defaultHelp(true)
				.version("${prog} " + VERSION)
				.description("DESCRIPTION\n"
+ "Soft clip reads on 5' and 3' ends.\n"
+ "BS-Seq libraries often have a bias at the 3' and/or 5' end of reads so it useful\n"
+ "to ignore these bases for methylation calling.\n"
+ "For source code see https://github.com/dariober/Java-cafe/tree/master/SoftClipBam"
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
					+ "Use - to print BAM to stdout or '-.sam' for sam to stdout.");
		
		parser.addArgument("--clipRead1", "-r1")
			.nargs(2)
			.type(Integer.class)
			.setDefault(0, 0)
			.help("A list of two integers: From read 1 soft clip this many bases from 5' and 3', respectively.");
		
		parser.addArgument("--clipRead2", "-r2")
			.nargs(2)
			.type(Integer.class)
			.setDefault(0, 0)
			.help("A list of two integers: From read 2 soft clip this many bases from 5' and 3', respectively.");

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
