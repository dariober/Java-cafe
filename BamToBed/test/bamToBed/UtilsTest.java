package bamToBed;

import static org.junit.Assert.*;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class UtilsTest {

	@Test
	public void canConvertSAMRecordToBed() {
		SAMRecord rec= new SAMRecord(null);
		rec.setReferenceName("chr1");
		rec.setAlignmentStart(1);
		rec.setCigarString("100M");
		rec.setReadName("read1");
		rec.setMappingQuality(30);
		assertEquals("chr1\t0\t100\tread1\t30\t+", Utils.SAMRecordToBed(rec));
	}

	@Test
	public void canFilterSAMRecord() {
		SAMRecord rec= new SAMRecord(null);
		rec.setMappingQuality(30);
		rec.setDuplicateReadFlag(true);
		rec.setReadNegativeStrandFlag(true);
		assertEquals(false, Utils.filterSamRecord(rec, 0, 0, 0));
		assertEquals(true, Utils.filterSamRecord(rec, 0, 0, 100));
		assertEquals(true, Utils.filterSamRecord(rec, 1025, 0, 0));
		assertEquals(true, Utils.filterSamRecord(rec, 0, 1024, 0));
		
	}
	
	
}

