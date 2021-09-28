package covariance.algorithms;

import java.io.*;
import java.util.*;

import covariance.datacontainers.*;
import covariance.utils.*;

public class JavaSCA implements ScoreGenerator
{
	private Alignment alignment;
	private double[][] locklessScores;
	private JavaSCA lastSubLockless;
	private int lastI = -99;
	char[] mostFrequentresidues;
	
	// from the Ranganathan lab
	public static float background[] =
	{0.072658f, 0.024692f, 0.050007f, 0.061087f,
        0.041774f, 0.071589f, 0.023392f, 0.052691f, 0.063923f,
        0.089093f, 0.023150f, 0.042931f, 0.052228f, 0.039871f,
        0.052012f, 0.073087f, 0.055606f, 0.063321f, 0.012720f,
        0.032955f}; 
	
	
	public JavaSCA( Alignment alignment ) throws Exception
	{
		this.alignment = alignment;
		initializeLocklessScores();
				
		mostFrequentresidues= Alignment.getMostFrequentResidues(alignment);
	}
	
	public double[][] getLocklessScores()
	{
		return this.locklessScores;
	}
	
	public double getLocklessSum(int i)
	{
		double sum = 0;
		
		for ( int x=0; x< MapResiduesToIndex.NUM_VALID_RESIDUES; x++ ) 
		{
			double temp = locklessScores[i][x];
			sum += temp * temp;
		}
		
		return Math.sqrt( sum );
	}
	
	private JavaSCA( Alignment fullAlignment, char[] residues, int i ) throws Exception
	{
		if ( residues[i] != '-' ) 
		{
			this.alignment = getSubsetAlignment( fullAlignment, i, residues[i] );
			initializeLocklessScores();
		}
	}
	
	public static Alignment getSubsetAlignment( Alignment fullAlignment, 
										   int filterResidueNumber, 
										   char filterResidueChar ) 
		throws Exception
	{
		if ( ! MapResiduesToIndex.isValidResidueChar( filterResidueChar ) ) 
			throw new Exception("Error!  " + filterResidueChar + " is not a vaild amino acid residue");
		
		if ( filterResidueNumber < 0 ) 
			throw new Exception("Error!  Filter residue number must be positive!");
		
		if ( filterResidueNumber > fullAlignment.getNumColumnsInAlignment() ) 
			throw new Exception("Error!  Filter residue number is longer than the sequence alignment");
		
		List subsetAlignmentLines = new ArrayList();
		
		for ( Iterator i = fullAlignment.getAlignmentLines().iterator();
			  i.hasNext(); ) 
		{
			AlignmentLine aLine = (AlignmentLine) i.next();
			
			if ( aLine.getSequence().charAt( filterResidueNumber ) == filterResidueChar ) 
				subsetAlignmentLines.add( aLine );
		}
		
		return new Alignment( filterResidueNumber + "_" + filterResidueChar + fullAlignment.getAligmentID(),
							 subsetAlignmentLines) ;
		
	}
	
	private static double getLnProbability(int aaIndex,
											  int N,
											  double nx) throws Exception
	{
		double returnVal = TTest.lnfactorial(N);
        returnVal -= TTest.lnfactorial(nx);
      	returnVal -= TTest.lnfactorial(N-nx);
      	returnVal += nx * Math.log(background[aaIndex]);
      	returnVal += (N-nx) * Math.log(1.0-background[aaIndex]);
      	
      	return returnVal;
	}
		
	private void initializeLocklessScores() throws Exception
	{
		int[][] counts = this.alignment.getCounts(); // initialize count array
		float[] frequencies = this.alignment.getFrequencies(); // initialize frequency array
		
		this.locklessScores = new double[ counts.length][MapResiduesToIndex.NUM_VALID_RESIDUES];
		
		for ( int x=0; x < counts.length; x++ ) 
		{	
			for ( int y=0; y< MapResiduesToIndex.NUM_VALID_RESIDUES; y++ ) 
			{
				this.locklessScores[x][y] = getLnProbability( y, 100, (((double)100*counts[x][y])/
							this.alignment.getNumSequencesInAlignment()));	
			}
		}
	}
	
	/**  This will work a lot faster if you call all the j's for the same i since the last
	 *   called subalignment i is cached
	 */
	public double getScore( Alignment alignment, int i, int j ) throws Exception
	{
		if ( mostFrequentresidues[i] == '-' ) 
			throw new Exception("Error!  No valid residues at position " + i );
		
		if ( this.alignment != alignment ) 
			throw new Exception( "Error!  Call with same alignment passed in from constructor");	
		
		
		if ( lastSubLockless == null || lastI != i  )
		{
			this.lastSubLockless = new JavaSCA( alignment, mostFrequentresidues, i );
			lastI = i;
		}
		
		double score = 0;
		
		for ( int x=0; x< MapResiduesToIndex.NUM_VALID_RESIDUES; x++ )
		{
			double temp = lastSubLockless.locklessScores[j][x] - this.locklessScores[j][x];
			score += temp * temp;
		}
		
		double returnVal = Math.sqrt( score );
		
		return returnVal;

	}
	
	public static void main(String[] args) throws Exception
	{
		if ( args.length != 2 ) 
		{
			System.out.println( "Usage JavaSCA inAlignment outFile" );
			return;
		}	
		
		Alignment a = new Alignment("1", new File( args[0]), false );
		
		JavaSCA sca = new JavaSCA(a);
		
		BufferedWriter writer = new BufferedWriter( new FileWriter( new File(
						args[1] )));
		
		writer.write( "i\tj\tscore\n");
		
		for ( int i =0; i < a.getNumColumnsInAlignment(); i++ ) 
			if ( a.columnHasValidResidue(i) ) 
				for ( int j = i + 1; j < a.getNumColumnsInAlignment(); j++ )
					if ( a.columnHasValidResidue(j) ) 
					writer.write( i + "\t" + j + "\t" + sca.getScore(a, i, j) + "\n" );				
		
		writer.flush();  writer.close();			
			
	}
	
	public boolean isSymmetrical()
	{
		return false;
	}
	
	public boolean reverseSort()
	{
		return false;
	}

	public String getAnalysisName()
	{
		return "JavaSCA";
	}
	
}
