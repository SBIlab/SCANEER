package covariance.datacontainers;

import covariance.parsers.PFamPdbAnnotationParser;

public class AlignmentPdbMatch implements Comparable
{
	private PFamPdbAnnotationParser parser;
	private float percentIdentity;
	private AlignmentLine aLine;
	private int shortestChainLength;
		
	public AlignmentPdbMatch( PFamPdbAnnotationParser parser, 
					 float percentIdentity, 
					 int shortestChainLength,
					 AlignmentLine aLine ) 
	{
		this.percentIdentity = percentIdentity;	
		this.aLine = aLine;
		this.parser = parser;
		this.shortestChainLength = shortestChainLength;
	}			
	
	
		
	public int compareTo(Object o)
	{
		if ( this.shortestChainLength!=  ((AlignmentPdbMatch)o).shortestChainLength) 
			return Float.compare( ((AlignmentPdbMatch)o).shortestChainLength, this.shortestChainLength);
				
			return Float.compare(  ((AlignmentPdbMatch)o).percentIdentity, this.percentIdentity);
	}

	public AlignmentLine getAlignmentLine()
	{
		return aLine;
	}

	public PFamPdbAnnotationParser getParser()
	{
		return parser;
	}

	public float getPercentIdentity()
	{
		return percentIdentity;
	}

	public int getShortestChainLength()
	{
		return shortestChainLength;
	}

	public String toString()
	{
		return aLine.getIdentifier() + " " + percentIdentity + "%  Length=" + shortestChainLength ;
	}

}
