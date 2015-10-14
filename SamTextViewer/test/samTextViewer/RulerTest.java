package samTextViewer;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

public class RulerTest {

	@Test
	public void testContructor() {
		
		int from= 1;
		int to= 10;
		int windowSize= 100;
		Ruler_TO_BE_DEPRECTED ruler= new Ruler_TO_BE_DEPRECTED(from, to, windowSize);
		assertEquals(1, ruler.getBpPerScreenColumn(), 0.001);

		ruler= new Ruler_TO_BE_DEPRECTED(1, 1000, 10);
		assertEquals(100, ruler.getBpPerScreenColumn(), 0.001);
		
		ruler= new Ruler_TO_BE_DEPRECTED(1, 160, 50);
		assertEquals(50, (int)ruler.getMapping().size());
	
		ruler= new Ruler_TO_BE_DEPRECTED(1, 100000, 170);
		assertEquals(170, (int)ruler.getMapping().size());
	}

	@Test
	public void canMapGenomePosToScreen() {
	
		Ruler_TO_BE_DEPRECTED ruler= new Ruler_TO_BE_DEPRECTED(100, 1000, 100);
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
		int windowSize= 333333;
		Ruler_TO_BE_DEPRECTED ruler= new Ruler_TO_BE_DEPRECTED(from, to, windowSize);
		// System.out.println(ruler.printableRuler(markDist));
		// assertEquals("50000000  50312480  50624960  50937440", ruler.printableRuler(markDist));
	}	
}
