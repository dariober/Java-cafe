package bamToBed;

import java.io.File;

import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class Utils {

	/**
	 * Open sam or bam file given filename or "-" to read from stdin.
	 * @param insam
	 * @param validationStringency See options in ValidationStringency. Typically one of:
	 * 		  ValidationStringency.SILENT, ValidationStringency.LENIENT, ValidationStringency.STRICT 
	 * @return
	 */
	public static SamReader reader(String insam, ValidationStringency validationStringency){
		SamReaderFactory sf = SamReaderFactory.
				makeDefault().
				validationStringency(validationStringency);
		SamReader sam= null;
		if(insam.equals("-")){
			SamInputResource resource= SamInputResource.of(System.in);
			sam= sf.open(resource);			
		} else {
			sam= sf.open(new File(insam));
		}
		return sam;
	}
	
}
