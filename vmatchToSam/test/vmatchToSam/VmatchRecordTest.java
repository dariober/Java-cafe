package vmatchToSam;

import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import org.junit.Test;

import exceptions.InvalidVmatchRecordException;
import htsjdk.samtools.SAMRecord;

public class VmatchRecordTest {

/*
# args=-showdesc 0 -s -p -d -complete -e 2 -q reads.fa /Users/berald01/svn_git/Java-cafe/trunk/vmatchToSam/test_data/ref.fa
25   ref_1   2   D 27   r1   0   2    7.44e-10  46    92.59
Sbjct: --TCTCTAGAGGTTAGCGAGACCGAAC                                           27
       !!                         
Query: CCTCTCTAGAGGTTAGCGAGACCGAAC                                           27


26   ref_1   1   D 27   r1   0   1    5.91e-12  50    96.30
Sbjct: -CTCTCTAGAGGTTAGCGAGACCGAAC                                           27
       !                          
Query: CCTCTCTAGAGGTTAGCGAGACCGAAC                                           27


27   ref_1   0   D 27   r1   0   0    1.84e-14  54   100.00
Sbjct: CCTCTCTAGAGGTTAGCGAGACCGAAC                                           27
Query: CCTCTCTAGAGGTTAGCGAGACCGAAC                                           27
*/
	
	@Test
	public void canConstructVmatchRecord() throws IOException, InvalidVmatchRecordException {
		BufferedReader br= new BufferedReader(new FileReader(new File("test_data/aln_1.vmatch.txt")));
		// First record
		VmatchRecord vRec= new VmatchRecord(br);
		assertEquals(25, vRec.getReferenceAlignmentLength());
		assertEquals(92.59, vRec.getPctIdentity(), 0.001);
		assertEquals("--TCTCTAGAGGTTAGCGAGACCGAAC", vRec.getAlignedReferenceSequence());
		assertEquals("CCTCTCTAGAGGTTAGCGAGACCGAAC", vRec.getAlignedQuerySequence());
		
		// Second record
		vRec= new VmatchRecord(br);
		assertEquals(26, vRec.getReferenceAlignmentLength());
		assertEquals(96.30, vRec.getPctIdentity(), 0.001);
		assertEquals("-CTCTCTAGAGGTTAGCGAGACCGAAC", vRec.getAlignedReferenceSequence());
		assertEquals("CCTCTCTAGAGGTTAGCGAGACCGAAC", vRec.getAlignedQuerySequence());
	}

	@Test
	public void canParseFile() throws IOException, InvalidVmatchRecordException {
		BufferedReader br= new BufferedReader(new FileReader(new File("test_data/aln_1.vmatch.txt")));
		VmatchRecord vRec= new VmatchRecord();
		int n= 0;
		while(br.ready()){
			vRec= new VmatchRecord(br);
			n++;
		}
		assertEquals(16, n);
	}
	
	@Test
	public void canGetSamRecordFromVmatchRecord() throws IOException, InvalidVmatchRecordException {
		BufferedReader br= new BufferedReader(new FileReader(new File("test_data/cigar.vmatch.txt")));
		VmatchRecord vRec= new VmatchRecord(br);
		SAMRecord samRec= vRec.getSAMRecord();
		assertEquals("4S20=2X22=4I32=", samRec.getCigar().toString());
		System.out.println(samRec.getSAMString());
	}	
}