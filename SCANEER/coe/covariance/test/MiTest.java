package covariance.test;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import covariance.algorithms.MICovariance;
import covariance.datacontainers.Alignment;
import covariance.datacontainers.AlignmentLine;
import covariance.utils.MapResiduesToIndex;

import junit.framework.TestCase;

public class MiTest extends TestCase
{

	public void testDistro1() throws Exception
	{
		List alignmentLines = new ArrayList();
		
		for ( int x=0; x < MapResiduesToIndex.NUM_VALID_RESIDUES; x++ ) 
			alignmentLines.add( new AlignmentLine("" + x, "" + 
													MapResiduesToIndex.getChararcter(x) + 
													MapResiduesToIndex.getChararcter(x) ));
				
				
		for ( int x=0; x< 10; x++)
			alignmentLines.add( new AlignmentLine("1", "-A"));
			
		for ( int x=0; x< 10; x++)
			alignmentLines.add( new AlignmentLine("1", "G-"));
		
													
		Alignment a = new Alignment("1", alignmentLines );
		MICovariance mi  = new MICovariance(a);
		
		
		double expectedScore = 20 * 0.05 *  Math.log(0.05/ (0.05*0.05)); 
		
		assertEquals( expectedScore, mi.getScore(a,0,1), 0.01);
		assertEquals( expectedScore, mi.getScore(a,1,0), 0.01);		
	}
	
	public void testOneConserved() throws Exception
	{
		List alignmentLines = new ArrayList();
		
		for ( int x=0; x < MapResiduesToIndex.NUM_VALID_RESIDUES; x++ ) 
			alignmentLines.add( new AlignmentLine("" + x, "Y" + 
													MapResiduesToIndex.getRandomChar() ));
		Alignment a = new Alignment("1", alignmentLines );
		MICovariance mi= new MICovariance(a);
		
		assertEquals( mi.getScore(a,0,1), 0, 0.01);	
	}
	
	public void testRandom() throws Exception
	{
		List alignmentLines = new ArrayList();
		
		for ( int x=0; x < MapResiduesToIndex.NUM_VALID_RESIDUES; x++ ) 
			alignmentLines.add( new AlignmentLine("" + x, "" + MapResiduesToIndex.getRandomChar() + 
													MapResiduesToIndex.getRandomChar() ));
		Alignment a = new Alignment("1", alignmentLines );
		MICovariance mi = new MICovariance(a);
		
		
		String iString = a.getColumnAsString(0);
		String jString= a.getColumnAsString(1);
		
		if ( iString.length() != jString.length() )
			throw new Exception("Logic error");
		
		StringBuffer iStringNew = new StringBuffer();
		StringBuffer jStringNew= new StringBuffer();
		
		for ( int x=0; x< iString.length(); x++ ) 
		{
			char iListChar = iString.charAt(x);;
			char jListChar = jString.charAt(x);
			
			if ( MapResiduesToIndex.isValidResidueChar( iListChar ) &&
				 MapResiduesToIndex.isValidResidueChar( jListChar ) )
			{
				iStringNew.append( iListChar );
				jStringNew.append( jListChar );
			}
		}
		
		iString = iStringNew.toString();
		jString = jStringNew.toString();
		int[] iFrequencies = mi.getFrequencies( iString );
		int[] jFrequencies = mi.getFrequencies( jString);
		
		
		for ( int x=0; x<MapResiduesToIndex.NUM_VALID_RESIDUES; x++ ) 
		{
			// no gaps so this should work
			assertEquals( iFrequencies[x],a.getCounts()[0][x]);
			assertEquals( jFrequencies[x],a.getCounts()[1][x]);
		}
		
		
		HashMap pairs = mi.getPairs( iString, jString, iFrequencies, jFrequencies );
		
		for ( Iterator i = pairs.keySet().iterator();
				i.hasNext(); )
		{
			String pairString = i.next().toString();
			
			MICovariance.PairsClass pClass = (MICovariance.PairsClass) pairs.get( pairString );
			
			assertEquals( pClass.pair, pairString );
			
			int count =0;
			
			for  (int x =0; x< a.getNumSequencesInAlignment(); x++ ) 
			{
				String alignmentSequence = ((AlignmentLine) a.getAlignmentLines().get(x)).getSequence();
				
				if ( alignmentSequence.charAt(0) == pairString.charAt(0) &&
						alignmentSequence.charAt(1) == pairString.charAt(1) ) 
						count++;
						
						
			}
			
			assertEquals( pClass.num, count );
		}
	}

	public void testAnAlignment() throws Exception
	{
		List alignmentLines = new ArrayList();
		
		alignmentLines.add( new AlignmentLine( "1",  "YH"));
		alignmentLines.add( new AlignmentLine( "1",  "YH"));
		alignmentLines.add( new AlignmentLine( "1",  "YH"));
		alignmentLines.add( new AlignmentLine( "1",  "YH"));
		alignmentLines.add( new AlignmentLine( "1",  "SF"));
		alignmentLines.add( new AlignmentLine( "1",  "SM"));
		alignmentLines.add( new AlignmentLine( "1",  "SM"));
		
		Alignment a = new Alignment("1", alignmentLines);	
		MICovariance mi = new MICovariance(a);
		
		assertEquals(mi.getScore(a, 0,1), 0.682908105, 0.00001);	
		assertEquals(mi.getScore(a, 1,0), 0.682908105, 0.00001);	
	}
	
	public void testAnotherAlignment() throws Exception
	{
		List alignmentLines = new ArrayList();
		
		alignmentLines.add( new AlignmentLine( "1",  "AC"));
		alignmentLines.add( new AlignmentLine( "1",  "AC"));
		alignmentLines.add( new AlignmentLine( "1",  "AC"));
		alignmentLines.add( new AlignmentLine( "1",  "AC"));
		alignmentLines.add( new AlignmentLine( "1",  "-C"));
		alignmentLines.add( new AlignmentLine( "1",  "AC"));
		alignmentLines.add( new AlignmentLine( "1",  "A-"));
		
		Alignment a = new Alignment("1", alignmentLines);	
		MICovariance mi = new MICovariance(a);
		
		assertEquals(mi.getScore(a, 0,1), 0.0, 0.001);	
		assertEquals(mi.getScore(a, 1,0), 0.0, 0.001);	
	}
	
	public void testYetAnotherOne() throws Exception
	{
		List alignmentLines = new ArrayList();
		
		alignmentLines.add( new AlignmentLine( "1",  "GCA"));
		alignmentLines.add( new AlignmentLine( "1",  "ACA"));
		alignmentLines.add( new AlignmentLine( "1",  "ACG"));
		alignmentLines.add( new AlignmentLine( "1",  "GGG"));
		alignmentLines.add( new AlignmentLine( "1",  "AGT"));
		alignmentLines.add( new AlignmentLine( "1",  "GTT"));
		alignmentLines.add( new AlignmentLine( "1",  "ATT"));
		alignmentLines.add( new AlignmentLine( "1",  "ATT"));
		alignmentLines.add( new AlignmentLine( "1",  "ATT"));
		alignmentLines.add( new AlignmentLine( "1",  "ATT"));
		
		Alignment a = new Alignment("1", alignmentLines);	
		MICovariance mi = new MICovariance(a);
		
		assertEquals(mi.getScore(a, 1,2),0.620686859, 0.001);	
		assertEquals(mi.getScore(a, 2,1), 0.620686859, 0.001);	
	}

	public MiTest(String arg0)
	{
		super(arg0);
	}

}
