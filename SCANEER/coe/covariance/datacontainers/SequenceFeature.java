package covariance.datacontainers;

public class SequenceFeature
{
	private int startInclusive;
	private int endInclusive;
	private String color;
	int sequenceNum;
	
	public SequenceFeature( int start, int end, String color, int sequenceNum ) 
	{
		this.startInclusive = start;
		this.endInclusive = end;
		this.color = color;
		this.sequenceNum = sequenceNum;
	}
	
	
	public String getColor()
	{
		return color;
	}

	public int getEndInclusive()
	{
		return endInclusive;
	}

	public int getStartInclusive()
	{
		return startInclusive;
	}

	public int getSequenceNum()
	{
		return sequenceNum;
	}

}
