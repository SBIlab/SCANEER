package covariance.test;


import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import covariance.algorithms.OmesCovariance;
import covariance.datacontainers.Alignment;
import covariance.datacontainers.AlignmentLine;
import covariance.datacontainers.AlignmentSubScore;
import covariance.utils.MapResiduesToIndex;

import junit.framework.TestCase;

public class OmesCovarianceTest extends TestCase
{
	public void testAnotherAlignment() throws Exception
	{
		List alignmentLines = new ArrayList();
		
		alignmentLines.add( new AlignmentLine( "1",  "YH"));
		alignmentLines.add( new AlignmentLine( "1",  "YH"));
		alignmentLines.add( new AlignmentLine( "1",  "YH"));
		alignmentLines.add( new AlignmentLine( "1",  "YH"));
		alignmentLines.add( new AlignmentLine( "1",  "TA"));
		alignmentLines.add( new AlignmentLine( "1",  "TA"));
		alignmentLines.add( new AlignmentLine( "1",  "TA"));
		alignmentLines.add( new AlignmentLine( "1",  "TA"));
		alignmentLines.add( new AlignmentLine( "1",  "TA"));
		alignmentLines.add( new AlignmentLine( "1",  "TL"));
		
		Alignment a = new Alignment("1", alignmentLines);	
		OmesCovariance oc= new OmesCovariance(a);
		assertEquals( oc.getScore(a,0,1), 1.98, 0.01);
		
		Collection subScores=  oc.getSubScores(a, 0,1);
		
		AlignmentSubScore subScore = AlignmentSubScore.getASubScore(subScores, 'Y', 'H');
		assertEquals(subScore.getNumExpected(), 1.6, 0.001 );
		assertEquals(subScore.getNumObserved(), 4);
		assertEquals(subScore.getScore(), 5.76, 0.01);
		
		subScore = AlignmentSubScore.getASubScore(subScores, 'T', 'L');
		assertEquals( subScore.getNumExpected(), 0.6, 0.01 );
		assertEquals( subScore.getNumObserved(), 1);
		assertEquals( subScore.getScore(), 0.16, 0.01 );
		
		
	}
	
	public void testYetAnotherAlignment() throws Exception
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
		OmesCovariance oc = new OmesCovariance(a);
		
		assertEquals(oc.getScore(a, 0,1), 1.306, 0.01);
		assertEquals(oc.getScore(a, 1,0), 1.306, 0.01);
		
		Collection subScores = oc.getSubScores(a, 0, 1);	
		
		AlignmentSubScore subScore = AlignmentSubScore.getASubScore(subScores, 'Y', 'H');
		assertEquals(subScore.getNumExpected(), 2.29, 0.01 );
		assertEquals(subScore.getNumObserved(), 4);
		assertEquals(subScore.getScore(), 2.94, 0.01);
		
		subScore = AlignmentSubScore.getASubScore(subScores, 'Y', 'F');
		assertEquals(subScore.getNumExpected(), 0.57, 0.01 );
		assertEquals(subScore.getNumObserved(), 0);
		assertEquals(subScore.getScore(), 0.33, 0.01);
		
		subScore = AlignmentSubScore.getASubScore(subScores, 'Y', 'M');
		assertEquals(subScore.getNumExpected(), 1.14, 0.01 );
		assertEquals(subScore.getNumObserved(), 0);
		assertEquals(subScore.getScore(), 1.31, 0.01);
		
		subScore = AlignmentSubScore.getASubScore(subScores, 'S', 'H');
		assertEquals(subScore.getNumExpected(), 1.71, 0.01 );
		assertEquals(subScore.getNumObserved(), 0);
		assertEquals(subScore.getScore(), 2.94, 0.01);
		
		subScore = AlignmentSubScore.getASubScore(subScores, 'S', 'F');
		assertEquals(subScore.getNumExpected(), 0.43, 0.01 );
		assertEquals(subScore.getNumObserved(), 1);
		assertEquals(subScore.getScore(), 0.33, 0.01);	
		
	}
	
	public void testKass() throws Exception
	{
		List alignmentLines = new ArrayList();
		
		alignmentLines.add( new AlignmentLine( "1",  "AHA"));
		alignmentLines.add( new AlignmentLine( "1",  "AHA"));		
		alignmentLines.add( new AlignmentLine( "1", "AHA"));
		alignmentLines.add( new AlignmentLine( "1",  "AHA"));
		alignmentLines.add( new AlignmentLine( "1",  "CFA"));
		alignmentLines.add( new AlignmentLine( "2",  "CGA"));
		alignmentLines.add( new AlignmentLine( "2",  "CGA"));
		alignmentLines.add( new AlignmentLine( "3",  "--A" ));
		alignmentLines.add( new AlignmentLine( "1",  "--A" ));
		alignmentLines.add( new AlignmentLine( "1",  "--A" ));
		alignmentLines.add( new AlignmentLine( "5",  "-L-" ));
		alignmentLines.add( new AlignmentLine( "6", "af-"));

	
		Alignment a = new Alignment("1", alignmentLines);	
		OmesCovariance oc = new OmesCovariance(a);
		assertEquals( oc.getScore(a,0,1), 1.31, 0.1);
		assertEquals( oc.getScore(a,0,2), oc.getScore(a, 2,0), 0.01);
		alignmentLines = new ArrayList();
		
		for ( int x=0; x< MapResiduesToIndex.NUM_VALID_RESIDUES; x++ ) 
		{
			alignmentLines.add( new AlignmentLine(""+x, 
						"" + MapResiduesToIndex.getChar(x) + MapResiduesToIndex.getChar(x)));
		}
		
		a = new Alignment("2", alignmentLines );
		assertEquals( new OmesCovariance(a).getScore(a,0,1), .95, 0.001);
	}

	public OmesCovarianceTest(String arg0)
	{
		super(arg0);
	}

}
