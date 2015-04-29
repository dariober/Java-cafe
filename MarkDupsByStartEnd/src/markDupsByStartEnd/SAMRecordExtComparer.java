package markDupsByStartEnd;

import java.util.Comparator;

public class SAMRecordExtComparer implements Comparator<SAMRecordExt> {
	
	@Override
	public int compare(SAMRecordExt x, SAMRecordExt y) {
	
		int startComparison = compare(x.getUnclippedStart(), y.getUnclippedStart());
		
		return startComparison != 0 ? startComparison
	                                : compare(x.getUnclippedEnd(), y.getUnclippedEnd());
	}
	
	private static int compare(int a, int b) {
		return a < b ? -1
				: a > b ? 1
				: 0;
	}
}
