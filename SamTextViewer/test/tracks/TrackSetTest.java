package tracks;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import exceptions.InvalidCommandLineException;
import samTextViewer.Utils;

public class TrackSetTest {

	@Test
	public void canSetYlimits() throws InvalidCommandLineException{
				
		String cmdInput= "#\\d+ 0 10 :ylim";
		// List<Track> tracks= new ArrayList<Track>();
		
		TrackSet ts= new TrackSet();
		
		Track t1= new Track(); t1.setFilename("foo.gz"); t1.setFileTag("#1"); ts.addOrReplace(t1);
		Track t2= new Track(); t2.setFilename("foo.txt"); t2.setFileTag("#20"); ts.addOrReplace(t2);
		Track t3= new Track(); t3.setFilename("bla.gz"); t3.setFileTag("#3"); ts.addOrReplace(t3);
		
		ts.setTrackYlimitsForRegex(cmdInput);
				
		assertEquals(0, ts.getTrackSet().get("#1").getYmin(), 0.001);
		assertEquals(0, ts.getTrackSet().get("#20").getYmin(), 0.001);
		assertEquals(10, ts.getTrackSet().get("#1").getYmax(), 0.001);
	}
	
}
