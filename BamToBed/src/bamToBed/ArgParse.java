package bamToBed;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.Arguments;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;

public class ArgParse {
	
	public static String PROG_NAME= "BamToBed";
	public static String VERSION= "0.1.0";
	
	/* Parse command line args */
	public static Namespace argParse(String[] args){
		ArgumentParser parser= ArgumentParsers
				.newArgumentParser(PROG_NAME)
				.defaultHelp(true)
				.version("${prog} " + VERSION)
				.description("DESCRIPTION\n"
+ "Convert BAM file to BED\n"
+ "See also https://github.com/dariober/Java-cafe/tree/master/BamToBed"
+ "\n"
+ "");	
		parser.addArgument("--inbam", "-i")
			.type(String.class)
			.required(false)
			.setDefault("-")
			.help("Input bam file. Should be sorted and indexed");

		parser.addArgument("--chrom", "-chrom")
			.type(String.class)
			.required(false)
			.setDefault("")
			.help("Only output reads on this chromosome");

		parser.addArgument("--from", "-from")
			.type(Integer.class)
			.required(false)
			.setDefault(0)
			.help("Only output reads starting from this position. 1-based, 0 for start of reference");
	
		parser.addArgument("--to", "-to")
			.type(Integer.class)
			.setDefault(0)
			.help("Only output reads up to this position. 1-based, 0 for end of reference");
		
		parser.addArgument("--mapq", "-q")
			.type(Integer.class)
			.required(false)
			.setDefault(0)
			.help("Minimum mapq score to output a read. Equivalent to samtools view -q");

		parser.addArgument("--requiredFlag", "-f")
			.type(Integer.class)
			.required(false)
			.setDefault(0)
			.help("Output reads with these bits set in flag. Equivalent to `samtools view -f`");

		parser.addArgument("--filterFlag", "-F")
			.type(Integer.class)
			.required(false)
			.setDefault(0)
			.help("Do not output reads with these bits set in flag. Equivalent to `samtools view -F`. Unmapped reads are always discared.");
		
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
