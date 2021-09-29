package covariance.algorithms;

public class RemainderClass implements Comparable
{
	
	private int orginalIndex;
	private float remainder;
	
	public int getOriginalIndex()
	{
		return this.orginalIndex;
	}
	
	public float getRemainder()
	{
		return this.remainder;
	}
	
	public RemainderClass( float remainder, int orginalIndex ) 
	{
		this.remainder = remainder;
		this.orginalIndex = orginalIndex;
	}
	
	/**  Sorts in ascending order
	 */
	public int compareTo(Object o)
	{
		return Float.compare( this.remainder, ((RemainderClass)o).remainder );
	}
	
	public String toString()
	{
		return orginalIndex + " " + remainder;
	}	
}
