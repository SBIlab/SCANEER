package covariance.algorithms;

import java.io.*;
import java.math.*;
import java.util.*;

import covariance.datacontainers.Alignment;
import covariance.utils.Factorials;
import covariance.utils.MapResiduesToIndex;

/**  Todo:  This and LocklessCoVariance could probably both become members of a Covariance subclass.
 *   ( Probably not really worth doing, though....)
 */
public class ELSCCovarianceSubAlignment
{
	private double[] dekkerCovarianceScores = null;
	private Alignment fullAlignment;
	private Alignment subsetAlignment;
	private static final BigInteger ONE = new BigInteger("1");
	
	public Alignment getFullAlignment()
	{
		return this.fullAlignment;
	}
	
	public Alignment getSubsetAlignment()
	{
		return this.subsetAlignment;
	}
	
	public List getSubsetAlignmentLines()
	{
		return subsetAlignment.getAlignmentLines();
	}
	
	public ELSCCovarianceSubAlignment(Alignment a, 
							  int filterResidueNum, 
							  char filterResidueChar) throws Exception
	{
		if ( filterResidueChar != '-' ) 
		{
			this.fullAlignment = a;
			this.subsetAlignment = JavaSCA.getSubsetAlignment( a, filterResidueNum, filterResidueChar );
			//System.out.println("Full = " + this.fullAlignment.getAlignmentLines().size() );
			//System.out.println("Subset = " + this.subsetAlignment.getAlignmentLines().size() );
		}
	}
	
	private int getSumAcrossResidues( int[][] counts, int position )
	{
		int sum = 0;
		
		for ( int x=0; x< MapResiduesToIndex.NUM_VALID_RESIDUES; x++ ) 
			sum+=counts[position][x];
	
		return sum;
	}
	
	
	/**  The sum of locklessIdealizedCounts should equal the total number of valid residues in that column
	 *   ( passed in as N )
	 */
	public static void adjustCounts( int[][] idealizedCounts, int N, LinkedList remainderList, int residuePosition )
		throws Exception
	{
		long numIdealized = 0;
		
		for ( int y=0; y < MapResiduesToIndex.NUM_VALID_RESIDUES; y++ ) 
			numIdealized+= idealizedCounts[residuePosition][y];
		
		if ( numIdealized == N ) 
			return;
		
		Collections.sort( remainderList );
		
		while ( numIdealized < N ) 
		{
			RemainderClass rClass = (RemainderClass) remainderList.removeLast();
			idealizedCounts[residuePosition][ rClass.getOriginalIndex()]++;
			numIdealized++;
		}
		
		// just a sanity check
		long check = 0;
		
		for ( int x=0; x< MapResiduesToIndex.NUM_VALID_RESIDUES; x++ ) 
			check+= idealizedCounts[residuePosition][x];
		
		if ( check != N )
			throw new Exception("Logic Error!");
	}
	
	
	private int[][] getIdealizedCounts(int[][] subsetCounts, int[][] fullCounts ) throws Exception
	{
		int[][] idealizedCounts = new int[ subsetCounts.length][MapResiduesToIndex.NUM_VALID_RESIDUES];
		
		for ( int x=0; x<subsetCounts.length; x++ ) 
		{
			LinkedList remainderList = new LinkedList();
			
			int n = getSumAcrossResidues( subsetCounts, x );
			int N = getSumAcrossResidues( fullCounts, x );
			
			for ( int y=0; y < MapResiduesToIndex.NUM_VALID_RESIDUES; y++ ) 
			{
				float idealCount = 0;
				
				if ( N > 0 ) 
					idealCount = ((float) n) * fullCounts[x][y] / N;
				
				idealizedCounts[x][y] = ( int ) idealCount;
				float remainder = idealCount - idealizedCounts[x][y];
				remainderList.add( new RemainderClass( remainder,y));
			}
			
			adjustCounts( idealizedCounts, n, remainderList, x);
		}
			
		return idealizedCounts;
	}
	
	public void dumpAllData( BufferedWriter writer ) throws Exception
	{
		double[] covarianceScores = getCovarianceScores();// initialize
		int[][] fullCounts = fullAlignment.getCounts();
		int[][] subsetCounts = subsetAlignment.getCounts();
		int[][] idealizedCounts = getIdealizedCounts( subsetCounts, fullCounts );
		
		writer.write("position\tresidue\tfullCounts\tsubsetCounts\tidealizedCounts\tcovarianceScores\n");
		
		for ( int x=0; x< fullAlignment.getNumColumnsInAlignment(); x++ ) 
			for ( int y=0; y < MapResiduesToIndex.NUM_VALID_RESIDUES; y++ ) 
			{
				writer.write( x + "\t" );
				writer.write( MapResiduesToIndex.getChar(y) + "\t" );
				writer.write( fullCounts[x][y] + "\t" );
				writer.write( subsetCounts[x][y] + "\t" );
				writer.write( idealizedCounts[x][y] + "\t" );
				writer.write( covarianceScores[x] + "\n" );
			}
		
		writer.flush();
		
	}
	
	public synchronized double[] getCovarianceScores() throws Exception
	{
		if ( dekkerCovarianceScores == null ) 
		{
			this.dekkerCovarianceScores = new double[fullAlignment.getNumColumnsInAlignment()];
			int[][] fullCounts = fullAlignment.getCounts();
			int[][] subsetCounts = subsetAlignment.getCounts();
			int[][] idealizedCounts = getIdealizedCounts( subsetCounts, fullCounts );
			
			for ( int x=0; x < fullAlignment.getNumColumnsInAlignment(); x++ ) 
			{
				double top = 0;
				double bottom = 0;
				
				for ( int y=0; y < MapResiduesToIndex.NUM_VALID_RESIDUES; y++ ) 
				{
					top += Factorials.getLnChoose( (int) fullCounts[x][y], (int) subsetCounts[x][y] );
						
					if ( fullCounts[x][y] > idealizedCounts[x][y]  ) 
							bottom += Factorials.getLnChoose( (int) fullCounts[x][y], (int)idealizedCounts[x][y]);
				}
				
				dekkerCovarianceScores[x] = - ( top - bottom);
			}
			
		}
		
		return dekkerCovarianceScores;
	}
}
