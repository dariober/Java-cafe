package tracks;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

public class TextProfileTest {

	@Test
	public void test() {
		
		List<Double> yValues= new ArrayList<Double>();
		yValues.add((double)0);
		yValues.add((double)0);
		yValues.add((double)0);
		yValues.add((double)0);
		
		TextProfile tp= new TextProfile(yValues, 10, Double.NaN, Double.NaN);
		System.out.println(tp.getProfile().size());
		System.out.println(tp.getProfile());
	}

}
