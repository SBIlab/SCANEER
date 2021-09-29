package covariance.datacontainers;

import java.util.Collection;
import java.util.Iterator;

public class PairColors
{
	private int xPos;
	private int yPos;
	private char xChar;
	private char yChar;
	private String colorString;
	
	public static String DEFAULT_COLOR_STRING= "black";
	
	public PairColors(int xPos, char xChar, int yPos, char yChar, String colorString)
	{
		this.xPos = xPos;
		this.xChar = xChar;
		this.yPos = yPos;
		this.yChar = yChar;
		this.colorString = colorString;	
	}
	
	public static String getColorStringAsRange( Collection pairColors, AlignmentLine aLine, int aPos )
	{
		String aSequence = aLine.getSequence();
		
		for ( Iterator i = pairColors.iterator();
				i.hasNext(); ) 
		{
			PairColors pColor = (PairColors) i.next();
			
			if ( aPos >= pColor.getXPos() && aPos <= pColor.getYPos()) 
			{
				return pColor.getColorString();
			}
		}
		
		return DEFAULT_COLOR_STRING;
		
	}
	
	public static String getColorString(  Collection pairColors, AlignmentLine aLine, int aPos )
	{
		String aSequence = aLine.getSequence();
		
		for ( Iterator i = pairColors.iterator();
				i.hasNext(); ) 
		{
			PairColors pColor = (PairColors) i.next();
			
			if ( aPos == pColor.getXPos() || aPos == pColor.getYPos()) 
			{
				if ( aSequence.charAt(pColor.getXPos()) == pColor.getXChar() &&
					 aSequence.charAt(pColor.getYPos() ) == pColor.getYChar() ) 
					 return pColor.getColorString();
			}
		}
		
		return DEFAULT_COLOR_STRING;
	
	}
	
	public boolean isMatch( int xPos, char xChar, int yPos, char yChar )	
	{
		if ( this.xPos != xPos ) 
			return false;
		
		if ( this.xChar != xChar )
			return false;
	
		if ( this.yPos != yPos ) 
			return false;
		
		if ( this.yChar != yChar ) 
			return false;
		
		return true;	
	}
	
	public static String getColorString( Collection pairColors, int xPos, char xChar, int yPos, char yChar )	
	{
		for ( Iterator i = pairColors.iterator();
				i.hasNext(); ) 
		{
			PairColors pColor = (PairColors) i.next();
			
			if ( pColor.isMatch(xPos, xChar, yPos, yChar))
				return pColor.getColorString();
		}
		
		return DEFAULT_COLOR_STRING;	
	}
	
	
	public String getColorString()
	{
		return colorString;
	}

	public char getXChar()
	{
		return xChar;
	}

	public int getXPos()
	{
		return xPos;
	}

	public char getYChar()
	{
		return yChar;
	}

	public int getYPos()
	{
		return yPos;
	}
}
