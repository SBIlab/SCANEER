package covariance.algorithms;

import covariance.datacontainers.*;


public interface ScoreGenerator
{
	public static final double NO_SCORE = -2;
	
	public double getScore( Alignment a, int i, int j ) throws Exception;
	
	public String getAnalysisName();
	
	/** return true is getScore(a,i,j) == getScore( a, j, i )
	 */
	public boolean isSymmetrical();
	
	/** return true if a low score indicates high covariance or conservation.
	 *  return false if a high score indicates high covariance or conservation.
	*/
	public boolean reverseSort();
	
}
