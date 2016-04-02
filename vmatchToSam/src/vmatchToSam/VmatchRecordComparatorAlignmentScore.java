package vmatchToSam;

import java.util.Comparator;

public class VmatchRecordComparatorAlignmentScore implements Comparator<VmatchRecord> {

	@Override
	public int compare(VmatchRecord v1, VmatchRecord v2) {
		
		// Descending order
		return v2.getAlignmentScore() - v1.getAlignmentScore();
		
	}

}
