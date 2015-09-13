package expandCigar;

import static org.junit.Assert.*;

import java.io.File;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import org.junit.Test;

import expandCigar.Utils;

public class UtilsTest {

	SamReaderFactory sf = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.STRICT);
	SamReader sam= sf.open(new File("test/test_data/header.sam"));
	SAMRecord rec= new SAMRecord(sam.getFileHeader());
	
	@Test
	public void testCanExpandMtoX() {

		rec.setReadBases("AAAAATCCCG".getBytes());
		rec.setCigarString("7M1I2M");
		rec.setAttribute("MD", "NN6N");
		String xcigar= Utils.expandCigarMtoX(rec);
		assertEquals("2X5=1I1=1X", xcigar.toString());
	}	
}
