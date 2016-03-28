package vmatchToSam;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.Arguments;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;

public class ArgParse {
	
	public static String PROG_NAME= "vmatchToSam";
	public static String VERSION= "0.1.0";
	
//	public static LinkedHashMap<String, String> docstrings= new LinkedHashMap<String, String>(); 
	
	
	/* Parse command line args */
	public static Namespace argParse(String[] args){
		ArgumentParser parser= ArgumentParsers
				.newArgumentParser(PROG_NAME)
				.defaultHelp(true)
				.version("${prog} " + VERSION)
				.description("DESCRIPTION\n"
+ "Convert output of vmatch to sam format.\n"
+ "For details see https://github.com/dariober/Java-cafe/tree/master/vmatchToSam\n"
+ "For vmatch aligner see http://www.vmatch.de/\n\n"
+ "Example\n"
+ "java vmatchToSam.jar vmatch.aln.txt\n"
+ "");	
		parser.addArgument("invmatch")
			.type(String.class)
			.required(false)
			.setDefault("-")
			.nargs("?")
			.help("Input vmatch files. Default or - to read from stdin");

		parser.addArgument("--fai", "-fai")
			.type(String.class)
			.required(false)
			.setDefault("")
			.help("Tab delimited file of reference sequence names and lengths. "
					+ "\n"
					+ "It can be the index of the reference sequences."
					+ "\n"
					+ "-fai or -fa must be supplied to produce the sam header.");
		
		parser.addArgument("--fasta", "-fa")
			.type(String.class)
			.required(false)
			.setDefault("")
			.help("Fasta file to extract sequence names and lengths.");

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
		
}
