package vmatchToSam;

import static org.junit.Assert.assertEquals;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.junit.Test;

import exceptions.InvalidVmatchRecordException;

public class VmatchRecordComparatorAlignmentScoreTest {

	@Test
	public void canSortRecordsByAlnScore() throws IOException, InvalidVmatchRecordException {
		BufferedReader br= new BufferedReader(new FileReader(new File("test_data/aln_1.vmatch.txt")));
		VmatchRecord vRec1= new VmatchRecord(br);
		VmatchRecord vRec2= new VmatchRecord(br);
		VmatchRecord vRec3= new VmatchRecord(br);
		
		vRec1.setAlignmentScore(10);
		vRec2.setAlignmentScore(20);
		vRec3.setAlignmentScore(15);
		
		List<VmatchRecord> recList= new ArrayList<VmatchRecord>(); 
		recList.add(vRec1);
		recList.add(vRec2);
		recList.add(vRec3);
		
		Collections.sort(recList, new VmatchRecordComparatorAlignmentScore());
		assertEquals(20, recList.get(0).getAlignmentScore());
		assertEquals(10, recList.get(2).getAlignmentScore());

	}

	
}
