package markDupsByStartEnd;

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
				.newArgumentParser("MarkDupsByStartEnd")
				.defaultHelp(true)
				.version("${prog} " + VERSION)
				.description("DESCRIPTION\n"
+ "Mark duplicate reads by looking at both the unclipped read start and end position."
+ "\n"
+ "For further details see"
+ "\n"
+ "https://github.com/dariober/Java-cafe/tree/master/MarkDupsByStartEnd"
+ "");	
		parser.addArgument("--insam", "-i")
			.type(String.class)
			.required(true)
			.help("Input sam or bam file. Use - to read from stdin.");
		
		parser.addArgument("--outsam", "-o")
			.type(String.class)
			.required(false)
			.setDefault("-")
			.help("Output file.\n"
					+ "Format will be sam or bam depending on extension.\n"
					+ "Use - to print SAM to stdout.");
		
		parser.addArgument("--unsortedOutput", "-us")
			.action(Arguments.storeTrue())
			.help("If set, output reads will be unsorted.\n"
					+ "By default reads are cooridnate sorted.");
		
		parser.addArgument("--ignoreReadGroup", "-rg")
			.action(Arguments.storeTrue())
			.help("Ignore read group info. If set, positional duplicates sharing *different* read groups "
					+ "will be considered duplicates.");
		
		parser.addArgument("--validationStringency", "-vs")
			.type(String.class)
			.required(false)
			.setDefault("SILENT")
			.choices("SILENT", "LENIENT", "STRICT")
			.help("Set picard validation stringency level.");
		
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
		
	}
}
