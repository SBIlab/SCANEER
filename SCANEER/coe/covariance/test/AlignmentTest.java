package covariance.test;

import java.util.ArrayList;
import java.util.List;

import covariance.datacontainers.Alignment;
import covariance.datacontainers.AlignmentLine;
import covariance.utils.MapResiduesToIndex;

import junit.framework.TestCase;

public class AlignmentTest extends TestCase
{

	public AlignmentTest(String arg0)
	{
		super(arg0);
	}
	
	
	public void testWeights() throws Exception
	{
		List list = new ArrayList();
		
		list.add( new AlignmentLine("0", "GYVGS"));
		list.add( new AlignmentLine("1", "GFDGF"));
		list.add( new AlignmentLine("2", "GYDGF"));
		list.add( new AlignmentLine("3", "GYQGG"));
	
		Alignment a= new Alignment( "1", list );
		
		double[] weights = a.adjustFrequenciesByWeights();
		
		assertEquals( weights.length, 4);
		assertEquals( weights[0], 0.267, 0.001);
		assertEquals( weights[1], 0.267, 0.001);
		assertEquals( weights[2], 0.200, 0.001);
		assertEquals( weights[3], 0.267, 0.001);
		
		char[] mostFrequentResidues = a.getMostFrequentResidues(a);
		
		assertEquals( mostFrequentResidues[0], 'G' );
		assertEquals( mostFrequentResidues[1], 'Y' );
		assertEquals( mostFrequentResidues[2], 'D');
		assertEquals( mostFrequentResidues[3], 'G' );
		assertEquals( mostFrequentResidues[4], 'F' );
	}
	
	public void testGetFilteredAlignment() throws Exception
	{
		List list = new ArrayList();
		
		list.add( new AlignmentLine("0", "AAAAAAAA"));
		list.add( new AlignmentLine("1", "GGGGGGGG"));
		list.add( new AlignmentLine("2", "HHHHHHHH"));
		list.add( new AlignmentLine("3", "AAAAAAAA"));
		list.add( new AlignmentLine("4", "---A----"));
		list.add( new AlignmentLine("5", "AAAAGGGG"));
		list.add( new AlignmentLine("6", "FFFFFAAA"));
		
		
		assertEquals( Alignment.getPairwiseIdentity( (AlignmentLine)list.get(0), (AlignmentLine) list.get(4)) , 
								100.0, 0.01 );	
								
		
		assertEquals( Alignment.getPairwiseIdentity( (AlignmentLine)list.get(0), (AlignmentLine) list.get(1)) , 
								0.0, 0.01 );	
								
		assertEquals( Alignment.getPairwiseIdentity( (AlignmentLine)list.get(0), (AlignmentLine) list.get(5)) , 
								50.0, 0.01 );	
								
		Alignment a = new Alignment("a", list );
		
		assertEquals(a.getNumSequencesInAlignment(), 7);
		
		char[] mostFrequentResidues = a.getMostFrequentResidues(a);
		
		assertEquals( mostFrequentResidues[0], 'A' );
		assertEquals( mostFrequentResidues[1], 'A' );
		assertEquals( mostFrequentResidues[6], 'A');
		
		a = a.getFilteredAlignment(50);
		
		assertEquals(a.getNumSequencesInAlignment(), 4);
		
		assertEquals( ((AlignmentLine)(a.getAlignmentLines().get(0))).getSequence(), "AAAAAAAA" );
		assertEquals( ((AlignmentLine)(a.getAlignmentLines().get(3))).getSequence(), "FFFFFAAA" );
		
		a = a.getFilteredAlignment(10);
		
		assertEquals(a.getNumSequencesInAlignment(), 3);
		assertEquals( ((AlignmentLine)(a.getAlignmentLines().get(0))).getSequence(), "AAAAAAAA" );
		assertEquals( ((AlignmentLine)(a.getAlignmentLines().get(2))).getSequence(), "HHHHHHHH" );
		
		/*  very slooooow!!
		a= PfamParser.getOneAlignment("ACT");
		
		a = a.getFilteredAlignment(20);
		
		for ( int x=0; x< a.getAlignmentLines().size(); x++)
		{
			AlignmentLine xLine = (AlignmentLine) a.getAlignmentLines().get(x);
			
			for ( int y=x+1; y < a.getAlignmentLines().size(); y++ ) 
			{
				AlignmentLine yLine = (AlignmentLine) a.getAlignmentLines().get(y);
				
				assertTrue( a.getPairwiseIdentity(xLine, yLine) < 20 );
			}
		} 
		
		a = a.getFilteredAlignment(0.00001f);
		
		assertEquals(a.getNumSequencesInAlignment(), 1);
		*/
	}

	public void testGetRatioValid() throws Exception
	{
		List list = new ArrayList();
		
		list.add( new AlignmentLine("1","AAA---G"));
		list.add( new AlignmentLine("2","AAAaaaG"));
		list.add( new AlignmentLine("3", "Cc----G"));
		list.add( new AlignmentLine("4", "F------"));
		
		Alignment a = new Alignment("", list);
		
		assertEquals( a.getRatioValid(0), 1.0, 0.1 );
		assertEquals( a.getRatioValid(1), .5, 0.1);
		assertEquals( a.getRatioValid(5), 0, 0.1);
		assertEquals( a.getRatioValid(6), 0, 0.75f);
	}

	public void testAlignmentConstructor() throws Exception
	{
		String someString = "AString";
		String someOtherString = "SomeOtherStringOFDifferentSize";
		
		List list = new ArrayList();
		
		list.add( new AlignmentLine("1", someString ));
		list.add( new AlignmentLine("2", someOtherString ));
		
		boolean threw = false;
		
		try
		{
			// should throw due to unequal alignment lengths
			Alignment a = new Alignment("1", list);
		}
		catch(Exception e) 
		{
			threw = true;
		}
		
		assertTrue( threw );
		
		String id = "The Name of this Alignment";
		list = TestUtils.getTestSequences();
		Alignment a = TestUtils.getTestAlignment( id, list );
		assertEquals( a.getAligmentID(), id);
		assertEquals( list.size(), a.getAlignmentLines().size());
		
		int length = list.iterator().next().toString().length();
		assertEquals( a.getNumColumnsInAlignment(), length);
	}

	public void testGetTotalNumValidResidues() throws Exception
	{
		String aString = TestUtils.getTestSequence();
		
		List list = new ArrayList();
		
		for ( int x=1; x <=20; x++ ) 
			list.add( new AlignmentLine( "" + x, aString ));
			
		StringBuffer nonValid = new StringBuffer();
		
		for ( int x=0; x< aString.length(); x++ ) 
			nonValid.append('-');
		
		String nonValidString = nonValid.toString();
			
		for ( int x=0; x< 14; x++ ) 
			list.add( new AlignmentLine( "" + (20+x), nonValidString ));
			
		Alignment a= new Alignment("1", list);
		
		int[] totalValid = a.getTotalNumValidResidues();
		
		assertEquals(totalValid.length, aString.length());
		
		for ( int x=0; x < MapResiduesToIndex.NUM_VALID_RESIDUES; x++ ) 
		{
			assertEquals(20 , totalValid[x]);			
		}
	}

	public void testGetCounts() throws Exception
	{
		StringBuffer buff = new StringBuffer();
		
		for ( int x=0; x < MapResiduesToIndex.NUM_VALID_RESIDUES; x++ ) 
			buff.append( "" + MapResiduesToIndex.getChararcter(x)  );
			
		String aSeq = buff.toString();
		
		List alignmentList = new ArrayList();		

		for ( int x=0; x <= 99; x++ ) 
			alignmentList.add( new AlignmentLine( ""+ x, aSeq + "-----" ));
			
		Alignment a= new Alignment("1", alignmentList);
		
		int count[][] = a.getCounts();
		
		for (int x=0; x< count.length; x++ ) 
			for ( int y=0; y< MapResiduesToIndex.NUM_VALID_RESIDUES; y++ ) 
			{
				if (  x==y ) 
					assertEquals(count[x][y], 100);
				else
					assertEquals(count[x][y], 0);
			}	
	}

	public void testGetFrequencies() throws Exception
	{
		List alignmentList = new ArrayList();		

		alignmentList.add( new AlignmentLine( "1", "----AAAAA---G" ));
		alignmentList.add( new AlignmentLine( "2", "----CCCCC---G" ));
		alignmentList.add( new AlignmentLine( "1", "----AAAAA---G" ));
		alignmentList.add( new AlignmentLine( "2", "----CCCCC----" ));
		alignmentList.add( new AlignmentLine( "1", "----AAAAA---G" ));
		alignmentList.add( new AlignmentLine( "2", "----CCCCC---G" ));
		alignmentList.add( new AlignmentLine( "1", "----AAAAA---G" ));
		alignmentList.add( new AlignmentLine( "2", "----CCCCC----" ));
			
		Alignment a= new Alignment("1", alignmentList);
		
		float[] frequencies = a.getFrequencies();
		
		for (int x=0; x< frequencies.length; x++ ) 
		{
			char c = MapResiduesToIndex.getChar(x);
			
			if ( c == 'A' || c == 'C'  ) 
				assertEquals( frequencies[x], (20f) / 46, 0.0001);
			else if ( c == 'G' ) 
				assertEquals(frequencies[x], 6f/46, 0.0001 );
			else
				assertEquals( frequencies[x], 0, 0.00001);	
		}	
		
		a = TestUtils.getTestAlignment("1", TestUtils.getTestSequences());
		
		float sum = 0;
		frequencies = a.getFrequencies();
		
		for ( int x=0; x < frequencies.length; x++ ) 
			sum += frequencies[x];
			
		assertEquals( 1, sum, 0.0005);
	}	

}
