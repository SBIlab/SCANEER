package covariance.datacontainers;

import covariance.utils.MapResiduesToIndex;

public class AlignmentLine
{
	public static final Character GAP_CHAR = new Character( '-');
	private String identifier;
	private String sequence;
	private int numValidChars;
	
	public String getIdentifier()
	{
		return identifier;
	}

	public String getSequence()
	{
		return sequence;
	}
	
	public AlignmentLine( String identifier, String sequence )
	{
		this.identifier = identifier;
		this.sequence = sequence;
		
		for ( int x=0; x< sequence.length(); x++ )
			if ( MapResiduesToIndex.isValidResidueChar( sequence.charAt(x) ) )
				numValidChars++;
	}	
	
	public String getUngappedSequence() throws Exception
	{
		StringBuffer buff = new StringBuffer();
		
		String capSequence = this.sequence.toUpperCase();
		
		for ( int x=0; x< capSequence.length(); x++ ) 
		{
			char c = capSequence.charAt( x );
			
			if ( MapResiduesToIndex.isValidResidueChar( c ) ) 
				buff.append(c);
		}
		
		return buff.toString();
	}
		
	public String toString()
	{
		return identifier + "\t" + sequence;
	}

	public int getNumValidChars()
	{
		return numValidChars;
	}

	public void setIdentifier(String identifier)
	{
		this.identifier = identifier;
	}

}
