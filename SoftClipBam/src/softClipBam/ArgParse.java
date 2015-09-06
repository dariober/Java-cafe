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
+ "Soft clip reads on 5'- and 3'-end.\n"
+ "BS-Seq libraries often have a bias at the 3' and/or 5' end of reads so it's\n"
+ "useful to ignore these bases for methylation calling.\n"
+ "\n"
+ "5' and 3' are based on the reference sequence. So for example these two reads\n"
+ "aligned to forward strand\n"
+ "HSX177:159:4:1104:10622:51324 99  chr18   10032   0  20M <-- Mate 1\n"
+ "HSX177:159:8:1107:24698:28804 163 chr18   10101   0  46M <-- Mate 2\n"
+ "\n"
+ "Clipped with\n"
+ "java -jar SoftClipBam.jar -i ear044_T3BS.chr18.bam -r1 0 10 -r2 20 0 \n"
+ "\n"
+ "Give:\n"
+ "HSX177:159:4:1104:10622:51324 99  chr18  10032  0  10M10S <-- Mate 1\n"
+ "HSX177:159:8:1107:24698:28804 163 chr18  10121  0  20S26M <-- Mate 2\n"
+ "\n"
+ "NB: Reads with no aligned bases, e.g. fully soft-clipped reads, are\n"
+ "returned as unmapped. This is consistent with the beahviour of htsjdk.\n"
+ "For source code see\n"
+ "https://github.com/dariober/Java-cafe/tree/master/SoftClipBam"
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
		
		parser.addArgument("--clipRead1", "-r1")
			.nargs(2)
			.type(Integer.class)
			.setDefault(0, 0)
			.help("A list of two integers: Soft clip read 1 these many\n"
					+ "bases from 5' and 3', respectively.");
		
		parser.addArgument("--clipRead2", "-r2")
			.nargs(2)
			.type(Integer.class)
			.setDefault(0, 0)
			.help("A list of two integers: Soft clip from read 2 these many\nbases from 5' and 3', respectively.");

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
