package covariance.test;

import java.util.ArrayList;
import java.util.List;

import covariance.datacontainers.Alignment;
import covariance.datacontainers.AlignmentLine;

import junit.framework.TestCase;

public class AlignmentFilterTest extends TestCase
{

	public void testColumnRemoval() throws Exception
	{
		List alignmentLines = new ArrayList();	
		
		alignmentLines.add( new AlignmentLine( "1", "AAC"));
		alignmentLines.add( new AlignmentLine( "2", "A-C"));
		alignmentLines.add( new AlignmentLine( "3", "AHY"));
		alignmentLines.add( new AlignmentLine( "4", "ADF"));
		
		Alignment a = new Alignment( "1", alignmentLines );
		boolean[] tossColumns = a.getTossColumns(0.7f);
		
		assertTrue( ! tossColumns[0] );
		assertTrue( tossColumns[1] );
		assertTrue( tossColumns[2] );
		
		tossColumns = a.getTossColumns(0.0f);
		
		assertTrue( ! tossColumns[0] );
		assertTrue( ! tossColumns[1] );
		assertTrue( ! tossColumns[2] );
		
		tossColumns = a.getTossColumns(0.4f);
		
		assertTrue( ! tossColumns[0] );
		assertTrue(  tossColumns[1] );
		assertTrue( ! tossColumns[2] );
		
		Alignment newA = a.getAlignmentWithRemovedUnconservedColumns(0.7f);
		
		List newLines = newA.getAlignmentLines();
		
		for ( int x=0; x<4; x++ ) 
		{
			AlignmentLine aLine = (AlignmentLine ) newLines.get(x);
			
			assertEquals(aLine.getSequence().length(), 1);
			assertEquals(aLine.getSequence().charAt(0), 'A');
		}
		
		newA = a.getAlignmentWithRemovedUnconservedColumns(0.45f);
		newLines = newA.getAlignmentLines();
		
		assertEquals( ((AlignmentLine) newLines.get(0)).getSequence(), "AC");
		assertEquals( ((AlignmentLine) newLines.get(1)).getSequence(), "AC");
		assertEquals( ((AlignmentLine) newLines.get(2)).getSequence(), "AY");
		assertEquals( ((AlignmentLine) newLines.get(3)).getSequence(), "AF");
		
		newA = a.getAlignmentWithRemovedUnconservedColumns(0.1f);
		newLines = newA.getAlignmentLines();
		
		assertEquals( ((AlignmentLine) newLines.get(0)).getSequence(), "AAC");
		assertEquals( ((AlignmentLine) newLines.get(1)).getSequence(), "A-C");
		assertEquals( ((AlignmentLine) newLines.get(2)).getSequence(), "AHY");
		assertEquals( ((AlignmentLine) newLines.get(3)).getSequence(), "ADF");
		
		
	}

	public AlignmentFilterTest(String arg0)
	{
		super(arg0);
	}

}
