package samTextViewer;

import static org.junit.Assert.*;

import org.junit.Test;

public class RulerTest {

	@Test
	public void testContructor() {
		
		int from= 1;
		int to= 10;
		int windowSize= 100;
		Ruler ruler= new Ruler(from, to, windowSize);
		assertEquals(1, ruler.bpPerScreenColumn(), 0.001);

		ruler= new Ruler(1, 1000, 10);
		assertEquals(100, ruler.bpPerScreenColumn(), 0.001);
		
		// Odd cases where expected window size differs from effective: 
		ruler= new Ruler(1, 160, 50);
		assertEquals(54, ruler.getScreenWidth());
		//System.out.println(ruler);
		ruler= new Ruler(1, 100000, 170);
		assertEquals(171, ruler.getScreenWidth());
	}

	@Test
	public void canMapGenomePosToScreen() {
	
		Ruler ruler= new Ruler(100, 1000, 100);
		// Corner cases:
		assertEquals(0, ruler.getScreenPositionAtGenomePosition(100));
		assertEquals(ruler.getMapping().size() - 1, ruler.getScreenPositionAtGenomePosition(1000));
		
		// First, in between, last:
		assertEquals(0, ruler.getScreenPositionAtGenomePosition(101));
		assertEquals(1, ruler.getScreenPositionAtGenomePosition(108));
		assertEquals(ruler.getMapping().size() - 1, ruler.getScreenPositionAtGenomePosition(999));
		
		// Out of range coords
		assertEquals(-1, ruler.getScreenPositionAtGenomePosition(99));
		assertEquals(-1, ruler.getScreenPositionAtGenomePosition(1001));
	}
	
	@Test
	public void canPrintRuler(){
		int from= 50000000;
		int to= 51000000;
		int markDist= 10;
		int windowSize= 33;
		Ruler ruler= new Ruler(from, to, windowSize);
		System.out.println(ruler.printableRuler(markDist));
		assertEquals("50000000  50294116  50588236  50882352", ruler.printableRuler(markDist));
	}	
}
