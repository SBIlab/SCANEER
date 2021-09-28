package covariance.test;

import java.math.BigDecimal;
import java.math.BigInteger;

import covariance.utils.Factorials;

import junit.framework.TestCase;

public class FactorialsTest extends TestCase
{

	public FactorialsTest(String arg0)
	{
		super(arg0);
	}
	

	public void testBinomialSum() throws Exception
	{
		assertEquals( Factorials.getBinomialSum(75, 0, 0.5f), 1, 0.01);
		
		BigDecimal bigD = Factorials.getBinomial(  75, 0, 0.5f);
		for ( int x=1; x<=30; x++ ) 
		{
			bigD = bigD.add(Factorials.getBinomial(  75, x, 0.5f));
		}
			
		assertEquals( bigD.doubleValue(), 1.0 - Factorials.getBinomialSum(75, 31, 0.5f), 0.01);
		
	}	

	public void testGetBinomial() throws Exception
	{
		BigDecimal aDistribution = Factorials.getBinomial(7,2,.583f);
		assertEquals( aDistribution.doubleValue(), 0.09, 0.01);
	}
	
	public void testLnChoose() throws Exception
	{
		
		assertEquals( Math.log(Factorials.choose(7,2).doubleValue()), 
				Factorials.getLnChoose(7,2) , 0.001);
				
		double lnTwoHundredFactorial = 859.6634733;
		double lnOneHundredNinetyEightFactorial = 849.0768721;
		
		double bottomLnchooseTwoHundredTwo = lnOneHundredNinetyEightFactorial + Math.log(2);
		double lnchooseTwoHundredTwo = lnTwoHundredFactorial - bottomLnchooseTwoHundredTwo;
		
		assertEquals(lnchooseTwoHundredTwo, Factorials.getLnChoose(200,2), 0.01);
		assertEquals(lnchooseTwoHundredTwo, Factorials.getLnChoose(200,198), 0.01);
		
		assertEquals( Factorials.getLnChoose(5000, 1000), Factorials.getLnChoose(5000, 4000), 0.01);
		
		assertEquals( Factorials.getLnBinomial(120, 23, 0.5f),
						Math.log(Factorials.getBinomial(120,23, 0.5f).doubleValue()), 0.01);
		

	}
	
	public void testFactorial() throws Exception
	{
		assertEquals( Factorials.getFactorial(5).intValue(), 120 );
		
		assertEquals( Factorials.choose( 5, 2).intValue(), 10);
		assertEquals( Factorials.choose( 6, 3).intValue(), 20 );

		assertEquals( Factorials.choose( 3,0).intValue(), 1);
		assertEquals( Factorials.choose( 0, 0).intValue(), 1);
		assertEquals( Factorials.choose( 1, 1).intValue(), 1);
		assertEquals( Factorials.choose( 30, 3).intValue(), Factorials.choose(30, 27).intValue() );
		
		assertEquals( Factorials.getFactorial(0).intValue(), 1);
		assertEquals( Factorials.getFactorial(1).intValue(), 1);
		assertEquals( Factorials.getFactorial(2).intValue(), 2);
		assertEquals( Factorials.getFactorial(3).intValue(), 6);
		
		for ( int x= Factorials.NUM_TO_CALCULATE-1; x > 0; x-- )
		{
			BigInteger newVal = Factorials.getFactorial(x).divide( Factorials.getFactorial(x-1));
			assertEquals( newVal.intValue(), x);
		}
	}


	public void testBinomial() throws Exception
	{
		BigDecimal sum = new BigDecimal("0");
		
		for ( int x=0; x<=150; x++) 
			sum = sum.add( Factorials.getBinomial(150, x, 0.5f));
		
		assertEquals(sum.doubleValue(), 1.0, 0.01);		
	}
}
