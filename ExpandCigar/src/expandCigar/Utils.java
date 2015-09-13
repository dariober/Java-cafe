package expandCigar;

import java.util.Arrays;
import java.util.List;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;

public class Utils {
	
	public static String expandCigarMtoX(SAMRecord rec) {
		byte[] refSeq= SequenceUtil.makeReferenceFromAlignment(rec, false);
		byte[] readSeq= rec.getReadBases();
		List<AlignmentBlock> alnBlocks= rec.getAlignmentBlocks();
		int i= 0;
		StringBuilder newCigarString= new StringBuilder(); 
		
		for(CigarElement cigar : rec.getCigar().getCigarElements()){
			if(cigar.getOperator().toString().equals("M")){
				AlignmentBlock block= alnBlocks.get(i);
				byte[] readBlock= Arrays.copyOfRange(readSeq, 
						block.getReadStart() - 1, 
						block.getReadStart() + block.getLength() - 1);
				byte[] refBlock= Arrays.copyOfRange(refSeq, 
						block.getReadStart() - 1, 
						block.getReadStart() + block.getLength() - 1);
				String ncigar= makeCigarStringFromReadAndRef(readBlock, refBlock);
				newCigarString.append(ncigar);
				i++;
			} else {
				Cigar xcigar= new Cigar();
				xcigar.add(cigar);
				newCigarString.append(xcigar);
			}
		}
		return newCigarString.toString();
	}

	private static String makeCigarStringFromReadAndRef(byte[] readSeq, byte[] refSeq){
		assert readSeq.length == refSeq.length;
		StringBuilder cigar = new StringBuilder();
		int len= 0;
		char current= '\0';
		char x= '\0';
		for(int i= 0; i < readSeq.length; i++){
			if(readSeq[i] == refSeq[i]){
				x= '=';
			} else {
				x= 'X';
			}
			if(current == '\0'){
				current= x;
			}
			if(current != x){
				cigar.append(len);
				cigar.append(current);
				len= 0;
				current= x;
			}
			len++;
		}
		cigar.append(len);
		cigar.append(current);
		return cigar.toString();
	}
	
	public static void print(Object x){
		String str= "";
		if (x instanceof byte[]){
			str= new String((byte[]) x);
		} else {
			str= x.toString();
		}
		System.out.println(str);
	}

}
