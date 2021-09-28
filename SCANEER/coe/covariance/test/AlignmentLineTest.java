package covariance.test;

import covariance.datacontainers.AlignmentLine;

import junit.framework.TestCase;

public class AlignmentLineTest extends TestCase
{

	public AlignmentLineTest(String arg0)
	{
		super(arg0);
	}

	public void testGetUngappedSequence() throws Exception
	{
	
		String testSequence = TestUtils.getTestSequence();	
		String prefix = "----El---vI-s------";
		String id = "someId";
		
		AlignmentLine aLine = new AlignmentLine( id, testSequence + prefix );
		assertEquals(aLine.getIdentifier(), id);
		assertEquals( aLine.getSequence(), testSequence + prefix );
		
		assertEquals( aLine.getUngappedSequence(), testSequence + "ELVIS" );
			
	}
	
	public void testNumValidChars() throws Exception
	{
		AlignmentLine aLine = new AlignmentLine("1","AAA---AA");
		assertEquals(aLine.getNumValidChars(), 5);
		
		aLine = new AlignmentLine("1", "abcdefg----------");
		assertEquals(aLine.getNumValidChars(), 0);
		
		aLine = new AlignmentLine("1", "FFFFFFFFFF" );
		assertEquals(aLine.getNumValidChars(), 10);
	}

}
