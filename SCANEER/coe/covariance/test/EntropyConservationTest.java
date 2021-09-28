package covariance.test;

import java.util.ArrayList;
import java.util.List;

import covariance.algorithms.ConservationSum;
import covariance.datacontainers.Alignment;
import covariance.datacontainers.AlignmentLine;
import covariance.utils.MapResiduesToIndex;

import junit.framework.TestCase;

public class EntropyConservationTest extends TestCase
{
	
	public void testEntropySomeMore() throws Exception
	{
		List list = new ArrayList();
		list.add( new AlignmentLine( "1", "AH" ));
		list.add( new AlignmentLine( "1", "AH" ));
		list.add( new AlignmentLine( "1", "CH" ));
		list.add( new AlignmentLine( "1", "CH" ));
		list.add( new AlignmentLine( "1", "GH" ));
		list.add( new AlignmentLine( "1", "GH" ));
		list.add( new AlignmentLine( "1", "YH" ));
		list.add( new AlignmentLine( "1", "YH" ));
		list.add( new AlignmentLine( "1", "YH" ));
				
		Alignment a = new Alignment("1", list);
		
		double score = (2f / 9f ) * Math.log(2f/9f);
		score = 3* score;
		
		score += ( 3f/ 9f ) * Math.log(3f/9f);
		score = -score;
		
		ConservationSum cSum = new ConservationSum(a);
		
		assertEquals(cSum.getScore(0), score, 0.01);
		assertEquals(cSum.getScore(1), 0, 0.001);
		assertEquals(cSum.getScore(a, 0,1), score /2, 0.01);
		assertEquals(cSum.getScore(a, 1,0), score /2, 0.01);	
	}
	
	public void testEntropy() throws Exception
	{
		List alignmentLines = new ArrayList();
		
		for ( int x=0; x<=10; x++ ) 
		{
			alignmentLines.add(new AlignmentLine("" + x, TestUtils.getTestSequence()));
		}
		
		Alignment a = new Alignment("1", alignmentLines );
		int length = TestUtils.getTestSequence().length();
		ConservationSum cSum = new ConservationSum(a);
		
		for ( int x=0; x< length; x++ )
			for ( int y=0; y < length; y++ ) 
			{
				assertEquals( cSum.getScore(a, x,y), 0.0, .0001);
				assertEquals( cSum.getScore(x), 0.0, .0001);
				assertEquals( cSum.getScore(y), 0.0, .0001);
			}
		
		alignmentLines = new ArrayList();
		
		for ( int x=0; x< MapResiduesToIndex.NUM_VALID_RESIDUES; x++ ) 
		{
			StringBuffer buff = new StringBuffer();
			
			for ( int y=0; y< MapResiduesToIndex.NUM_VALID_RESIDUES; y++ ) 
				buff.append( MapResiduesToIndex.getChar(x) );
				
			alignmentLines.add( new AlignmentLine("" + x, buff.toString()));
		}
		
		alignmentLines.add( new AlignmentLine("1", "aaaaaaaaaaaaaaaaaaaa"));
		alignmentLines.add( new AlignmentLine("1", "--------------------"));
		
		a = new Alignment("2", alignmentLines);
		cSum = new ConservationSum(a);
		double expectedScore = -20 * .05 * Math.log(.05);
		
		for ( int x=0; x < MapResiduesToIndex.NUM_VALID_RESIDUES; x++ )
			for ( int y=0; y < MapResiduesToIndex.NUM_VALID_RESIDUES; y++ ) 
			{
				if ( x != y ) 
					assertEquals( cSum.getScore(a, x,y), expectedScore, 0.01  );
			}
	}

	public EntropyConservationTest(String arg0)
	{
		super(arg0);
	}

}
