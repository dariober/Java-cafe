package samTextViewer;

import static org.junit.Assert.*;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

import org.junit.Test;

public class TextReadTest {

	@Test
	public void getNoAlignedBasesFromCigar() {
		// A reminder of how cigar works.
		SAMRecord rec= new SAMRecord(null);
		rec.setCigarString("2M2D4M");
		int noRefBases= rec.getCigar().getReferenceLength();
		assertEquals(8, noRefBases);
		
		rec.setCigarString("2M2N4M"); // N same as D
		noRefBases= rec.getCigar().getReferenceLength();
		assertEquals(8, noRefBases);
		
		rec.setCigarString("2M2I4M");
		noRefBases= rec.getCigar().getReferenceLength();
		assertEquals(6, noRefBases);
		
		rec.setCigarString("2S10M4H");
		noRefBases= rec.getCigar().getReferenceLength();
		assertEquals(10, noRefBases);		
		
		// Cigar oprators consuming:
		rec.setCigarString("2H"); // Hard clipping doesn't consume bases on read or ref
		CigarElement el= rec.getCigar().getCigarElement(0);
		assertFalse(el.getOperator().consumesReadBases());		
		assertFalse(el.getOperator().consumesReferenceBases());

		rec.setCigarString("2S"); // Clips the read not the ref.
		el= rec.getCigar().getCigarElement(0);
		assertTrue(el.getOperator().consumesReadBases());		
		assertFalse(el.getOperator().consumesReferenceBases());

		rec.setCigarString("2D"); // This makes a gap in the read when aligned to ref
		el= rec.getCigar().getCigarElement(0);
		assertFalse(el.getOperator().consumesReadBases());		
		assertTrue(el.getOperator().consumesReferenceBases());

		rec.setCigarString("2I"); // Opposite of D
		el= rec.getCigar().getCigarElement(0);
		assertTrue(el.getOperator().consumesReadBases());		
		assertFalse(el.getOperator().consumesReferenceBases());
		
	}

	@Test
	public void canGetReadStartEndInWindowCoords(){		
		/* Window
		100																			   180
		|                                                                              |
		AACCGGTT
		 */
		
		TextWindow textWindow= new TextWindow();
		textWindow.setFrom(100);
		textWindow.setTo(180);
		
		SAMRecord rec= new SAMRecord(null);
		rec.setAlignmentStart(100);
		rec.setCigarString("100M");
		TextRead textRead= new TextRead(rec, textWindow);
		assertEquals(1, textRead.getTextStart());
		
		rec.setAlignmentStart(90);
		textRead= new TextRead(rec, textWindow);
		assertEquals(1, textRead.getTextStart());

		rec.setAlignmentStart(101);
		textRead= new TextRead(rec, textWindow);
		assertEquals(2, textRead.getTextStart());

		// End of read beyond end of window.
		rec.setAlignmentStart(100);
		rec.setCigarString("100M");
		textRead= new TextRead(rec, textWindow);
		assertEquals(80, textRead.getTextEnd());

		rec.setAlignmentStart(101);
		rec.setCigarString("10M");
		textRead= new TextRead(rec, textWindow);
		assertEquals(10, textRead.getTextEnd());
	}
	
	@Test
	public void canGetTextReadInWindow(){

		TextWindow textWindow= new TextWindow();
		textWindow.setFrom(100);
		textWindow.setTo(180);
		
		SAMRecord rec= new SAMRecord(null);
		rec.setAlignmentStart(98);
		rec.setCigarString("2M2=2X");
		rec.setReadString("AACCGG");
		TextRead textRead= new TextRead(rec, textWindow);
		assertEquals("CCGG", new String(textRead.getTextRead()));

		rec= new SAMRecord(null);
		rec.setAlignmentStart(120);
		rec.setCigarString("6M");
		rec.setReadString("AACCGG");
		textRead= new TextRead(rec, textWindow);
		assertEquals("AACCGG", new String(textRead.getTextRead()));

		rec= new SAMRecord(null);
		rec.setAlignmentStart(180);
		rec.setCigarString("6M");
		rec.setReadString("AACCGG");
		textRead= new TextRead(rec, textWindow);
		assertEquals("A", new String(textRead.getTextRead()));

		rec= new SAMRecord(null);
		rec.setAlignmentStart(100);
		rec.setCigarString("3M4D4N3M");
		rec.setReadString("AACCGG");		
		textRead= new TextRead(rec, textWindow);
		assertEquals("AAC----____CGG", new String(textRead.getTextRead()));

		rec= new SAMRecord(null);
		rec.setAlignmentStart(100);
		rec.setCigarString("2M2I2M2S4H");
		rec.setReadString("AACCGGNN");		
		textRead= new TextRead(rec, textWindow);
		assertEquals("AAGG", new String(textRead.getTextRead()));		
	}
	
	@Test
	public void canGetPrintableRead(){

		TextWindow textWindow= new TextWindow();
		textWindow.setFrom(100);
		textWindow.setTo(110);

		SAMRecord rec= new SAMRecord(null);
		rec.setAlignmentStart(101);
		rec.setCigarString("6M");
		rec.setReadString("AACCGG");		
		TextRead textRead= new TextRead(rec, textWindow);
		assertEquals(" AACCGG", new String(textRead.getPrintableText()));

		textWindow.setTo(103);
		rec.setAlignmentStart(103);
		rec.setCigarString("6M");
		rec.setReadString("AACCGG");		
		textRead= new TextRead(rec, textWindow);
		assertEquals("   A", new String(textRead.getPrintableText()));

	}
	
	@Test
	public void canEncodeReadToMatchRef(){

		TextWindow textWindow= new TextWindow();
		textWindow.setFrom(100);
		textWindow.setTo(110);

		SAMRecord rec= new SAMRecord(null);
		rec.setReadNegativeStrandFlag(true);
		rec.setAlignmentStart(100);
		rec.setCigarString("6M3D2M");
		rec.setReadString("AACCGGTT");
		byte[] refSeq=    "AACCGNNNNNN".getBytes();		
		TextRead textRead= new TextRead(rec, textWindow);
		assertEquals(",,,,,g---tt", new String(textRead.getTextChars(refSeq)));

	}
	@Test
	public void canConvertReadToBS(){

		TextWindow textWindow= new TextWindow();
		textWindow.setFrom(100);
		textWindow.setTo(110);
		
		SAMRecord rec= new SAMRecord(null);
		byte[] refSeq=    "AACCGG".getBytes();
		rec.setReadString("AACTGA");
		rec.setAlignmentStart(100);
		rec.setCigarString("6M");
		TextRead textRead= new TextRead(rec, textWindow);
		
		// Unpaired reads 
		// --------------
		rec.setReadPairedFlag(false);
		
		rec.setReadNegativeStrandFlag(false);
		assertEquals("..MU.A", new String(textRead.convertTextToBS(refSeq)));

		rec.setReadNegativeStrandFlag(true);
		assertEquals(",,,tmu", new String(textRead.convertTextToBS(refSeq)));
		
		// Paired: First in pair (same as unpaired)
		// --------------------------------
		rec.setReadPairedFlag(true);
		rec.setSecondOfPairFlag(false);

		rec.setReadNegativeStrandFlag(false);
		assertEquals("..MU.A", new String(textRead.convertTextToBS(refSeq)));

		rec.setReadNegativeStrandFlag(true);
		assertEquals(",,,tmu", new String(textRead.convertTextToBS(refSeq)));		

		// Paired: Second in pair
		// ------------
		rec.setReadPairedFlag(true);
		rec.setSecondOfPairFlag(true);

		rec.setReadNegativeStrandFlag(false);	        
		assertEquals("...TMU", new String(textRead.convertTextToBS(refSeq)));
		
		rec.setReadNegativeStrandFlag(true);
		assertEquals(",,mu,a", new String(textRead.convertTextToBS(refSeq)));
	}
}
