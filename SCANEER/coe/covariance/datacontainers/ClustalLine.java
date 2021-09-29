package covariance.datacontainers;

public class ClustalLine
{
	private int positionInOriginatingSequence;
	private int clustalIndex;
	private boolean isAlignmentSequence;
	private String lineName;
		
		
	public ClustalLine( 
					 String lineName,
		 			 int clustalIndex, 
					 int positionInOriginatingSequence, 
					 boolean isAlignmentSequence ) 
	{
		this.lineName = lineName;
		this.clustalIndex = clustalIndex;
		this.positionInOriginatingSequence = positionInOriginatingSequence;
		this.isAlignmentSequence = isAlignmentSequence;
	}		
		
	public int getClustalIndex()
	{
		return clustalIndex;
	}

	public boolean isAlignmentSequence()
	{
		return isAlignmentSequence;
	}

	public String getLineName()
	{
		return lineName;
	}

	public int getPositionInOriginatingSequence()
	{
		return positionInOriginatingSequence;
	}

	public void incrementPositionInOriginalSequence()
	{
		this.positionInOriginatingSequence++;
	}
}
