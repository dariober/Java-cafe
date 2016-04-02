package utils;

import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;

import org.junit.Test;

import exceptions.InvalidVmatchRecordException;
import htsjdk.samtools.SAMFileHeader;
import vmatchToSam.VmatchRecord;

public class UtilsTest {

	@Test
	public void canCleanRefSeqName(){
		assertEquals("foo", Utils.refSeqNameForSam("  foo"));
		assertEquals("foo", Utils.refSeqNameForSam("foo bar"));
	}
	
	@Test
	public void canMakeHeaderFromFastaIndex() throws IOException {

		SAMFileHeader fh = Utils.makeSAMHeaderFromFastaIndex("test_data/ref.fa.fai");
		assertEquals("ref_1", fh.getSequence(0).getSequenceName());
		
	}

	@Test
	public void canMakeHeaderFromFasta() throws IOException {

		SAMFileHeader fh = Utils.makeSAMHeaderFromFasta("test_data/ref_noindex.fa");
		assertEquals("ref_1", fh.getSequence(0).getSequenceName());
		assertEquals(373, fh.getSequence(0).getSequenceLength());
		assertEquals("ref_2", fh.getSequence(1).getSequenceName());
		assertEquals(215, fh.getSequence(1).getSequenceLength());
	
	}	
	
}
