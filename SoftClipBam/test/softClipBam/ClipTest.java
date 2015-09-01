package softClipBamReads;

import static org.junit.Assert.*;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import java.io.File;
import org.junit.Test;

/* 
 
5' +++++>>>>R1>>>>+++++++++++++++++++++++++++++>>>>R2>>>>++++++++++++++++++++ 3'
3' ---------------<<<<R2<<<<--------------------------------<<<<R1<<<<------- 5'

*/

public class ClipTest {

	SamReaderFactory sf = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.STRICT);
	SamReader sam= sf.open(new File("test/test_data/header.sam"));

	@Test
	public void readUnchanged(){

		// INPUT R1 Single end forward
		SAMRecord rec= new SAMRecord(sam.getFileHeader());
		rec.setReferenceName("chr18");
		rec.setAlignmentStart(1000);
		rec.setReadBases("AAAAACCCCCTTTTTGGGGG".getBytes());		
		rec.setCigarString("20M");
		
		Clipper.clip(rec, 0, 0, 0, 0);
		assertEquals("20M", rec.getCigarString());

		// Single end reverse
		rec.setCigarString("20M");
		rec.setReadNegativeStrandFlag(true);
		Clipper.clip(rec, -1, -1, -1, -1);
		assertEquals("20M", rec.getCigarString());
		
		// Read 2 on reverse
		rec.setReadPairedFlag(true);
		rec.setSecondOfPairFlag(true);
		rec.setReadNegativeStrandFlag(true);
		Clipper.clip(rec, -1, -1, -1, -1);
		assertEquals("20M", rec.getCigarString());

		// Read 1 on forward
		rec.setReadPairedFlag(true);
		rec.setFirstOfPairFlag(true);
		rec.setReadNegativeStrandFlag(false);
		Clipper.clip(rec, -1, -1, -1, -1);
		assertEquals("20M", rec.getCigarString());
		
		// Read 1 on reverse
		rec.setReadPairedFlag(true);
		rec.setFirstOfPairFlag(true);
		rec.setReadNegativeStrandFlag(true);
		Clipper.clip(rec, -1, -1, -1, -1);
		assertEquals("20M", rec.getCigarString());		
	}

	@Test
	public void canClip3PrimeEnd() {
		
		// INPUT Single end forward
		SAMRecord rec= new SAMRecord(sam.getFileHeader());
		rec.setReadName("myread");
		rec.setReferenceName("chr18");
		rec.setAlignmentStart(1000);
		rec.setReadBases("AAAAACCCCCTTTTTGGGGG".getBytes());		
		
		// Single end forward
		rec.setCigarString("20M");
		Clipper.clip(rec, 0, 2, 0, 0);
		assertEquals("18M2S", rec.getCigarString());
		assertEquals(1000, rec.getAlignmentStart());
		
		// Single end reverse
		rec.setCigarString("20M");
		rec.setReadNegativeStrandFlag(true);
		Clipper.clip(rec, 0, 2, 0, 0);
		assertEquals("2S18M", rec.getCigarString());
		assertEquals(1000+2, rec.getAlignmentStart());

		// First in pair forward (same as SE)
		rec.setAlignmentStart(1000);
		rec.setCigarString("20M");
		rec.setReadPairedFlag(true);
		rec.setFirstOfPairFlag(true);
		rec.setReadNegativeStrandFlag(false);
		Clipper.clip(rec, 0, 2, 0, 0);
		assertEquals("18M2S", rec.getCigarString());
		assertEquals(1000, rec.getAlignmentStart());
		
		// First in pair reverse
		rec.setAlignmentStart(1000);
		rec.setCigarString("20M");
		rec.setReadPairedFlag(true);
		rec.setFirstOfPairFlag(true);
		rec.setReadNegativeStrandFlag(true);
		Clipper.clip(rec, 0, 2, 0, 0);
		assertEquals("2S18M", rec.getCigarString());
		assertEquals(1000+2, rec.getAlignmentStart());

		// ---------
		// Second in pair forward
		rec= new SAMRecord(sam.getFileHeader());
		rec.setReadName("myread");
		rec.setReferenceName("chr18");
		rec.setAlignmentStart(1000);
		rec.setReadBases("AAAAACCCCCTTTTTGGGGG".getBytes());		
		rec.setReadPairedFlag(true);
		rec.setSecondOfPairFlag(true);
		rec.setReadNegativeStrandFlag(false);

		rec.setCigarString("20M");
		Clipper.clip(rec, 0, 0, 0, 2);
		assertEquals("18M2S", rec.getCigarString());
		assertEquals(1000, rec.getAlignmentStart());

		// Second in pair reverse
		rec.setReadNegativeStrandFlag(true);
		rec.setCigarString("20M");
		Clipper.clip(rec, 0, 0, 0, 2);
		assertEquals("2S18M", rec.getCigarString());
		assertEquals(1000+2, rec.getAlignmentStart());
			}
	
	@Test
	public void canClip5PrimeEnd() {
		
		// INPUT Single end forward
		SAMRecord rec= new SAMRecord(sam.getFileHeader());
		rec.setReadName("myread");
		rec.setReferenceName("chr18");
		rec.setAlignmentStart(1000);
		rec.setReadBases("AAAAACCCCCTTTTTGGGGG".getBytes());		
		
		// Single end forward
		rec.setCigarString("20M");
		Clipper.clip(rec, 2, 0, 0, 0);
		assertEquals("2S18M", rec.getCigarString());
		assertEquals(1000+2, rec.getAlignmentStart());
		
		// Single end reverse
		rec.setAlignmentStart(1000);
		rec.setCigarString("20M");
		rec.setReadNegativeStrandFlag(true);
		Clipper.clip(rec, 2, 0, 0, 0);
		assertEquals("18M2S", rec.getCigarString());
		assertEquals(1000, rec.getAlignmentStart());

		// First in pair forward (same as SE)
		rec.setAlignmentStart(1000);
		rec.setCigarString("20M");
		rec.setReadPairedFlag(true);
		rec.setFirstOfPairFlag(true);
		rec.setReadNegativeStrandFlag(false);
		Clipper.clip(rec, 2, 0, 0, 0);
		assertEquals("2S18M", rec.getCigarString());
		assertEquals(1000+2, rec.getAlignmentStart());
		
		// First in pair reverse
		rec.setAlignmentStart(1000);
		rec.setCigarString("20M");
		rec.setReadPairedFlag(true);
		rec.setFirstOfPairFlag(true);
		rec.setReadNegativeStrandFlag(true);
		Clipper.clip(rec, 2, 0, 0, 0);
		assertEquals("18M2S", rec.getCigarString());
		assertEquals(1000, rec.getAlignmentStart());

		// ---------
		// Second in pair forward
		rec= new SAMRecord(sam.getFileHeader());
		rec.setReadName("myread");
		rec.setReferenceName("chr18");
		rec.setAlignmentStart(1000);
		rec.setReadBases("AAAAACCCCCTTTTTGGGGG".getBytes());		
		rec.setReadPairedFlag(true);
		rec.setSecondOfPairFlag(true);
		rec.setReadNegativeStrandFlag(false);

		rec.setCigarString("20M");
		Clipper.clip(rec, 0, 0, 2, 0);
		assertEquals("2S18M", rec.getCigarString());
		assertEquals(1000+2, rec.getAlignmentStart());

		// Second in pair reverse
		rec= new SAMRecord(sam.getFileHeader());
		rec.setReadName("myread");
		rec.setReferenceName("chr18");
		rec.setAlignmentStart(1000);
		rec.setReadBases("AAAAACCCCCTTTTTGGGGG".getBytes());		
		rec.setReadPairedFlag(true);
		rec.setSecondOfPairFlag(true);
		rec.setReadNegativeStrandFlag(true);

		rec.setCigarString("20M");
		Clipper.clip(rec, 0, 0, 2, 0);
		assertEquals("18M2S", rec.getCigarString());
		assertEquals(1000, rec.getAlignmentStart());
	
	}
	
	@Test
	public void canHandleSoftClippedReads() {

		SAMRecord rec= new SAMRecord(sam.getFileHeader());
		rec.setReadName("myread");
		rec.setReferenceName("chr18");
		rec.setAlignmentStart(1000);
		rec.setReadBases("AAAAACCCCCTTTTTGGGGG".getBytes());		
		rec.setCigarString("2S18M");
		Clipper.clip(rec, 2, 0, 0, 0);
		assertEquals("2S18M", rec.getCigarString());
		assertEquals(1000, rec.getAlignmentStart());
		
		rec.setCigarString("18M2S");
		Clipper.clip(rec, 0, 2, 0, 0);
		assertEquals("18M2S", rec.getCigarString());
		assertEquals(1000, rec.getAlignmentStart());
	}

	@Test
	public void canHandleIndels() {

		SAMRecord rec= new SAMRecord(sam.getFileHeader());
		rec.setReadName("myread");
		rec.setReferenceName("chr18");
		rec.setAlignmentStart(1000);
		rec.setReadBases("AAAAACCCCCTTTTTGGGGG".getBytes());		

		rec.setCigarString("1M2I17M");
		Clipper.clip(rec, 2, 0, 0, 0);
		assertEquals("3S17M", rec.getCigarString());	
		assertEquals(rec.getReadLength(), rec.getCigar().getReadLength());
		
	}
}
