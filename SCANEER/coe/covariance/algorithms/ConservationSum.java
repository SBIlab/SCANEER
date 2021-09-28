package covariance.algorithms;

import java.io.*;

import covariance.datacontainers.*;
import covariance.utils.*;

public class ConservationSum implements ScoreGenerator, ConservationGenerator
{
	private double[] absoluteSequenceEntropy;
	private Alignment alignment;
	private String id;
	
	public double getScore( Alignment a, int i, int j ) throws Exception
	{
		if ( a != this.alignment ) 
			throw new Exception("Please call getScore() with the same alignment used in the constructor");
		
		return (absoluteSequenceEntropy[i] + absoluteSequenceEntropy[j]) / 2;
	}
	
	public static void main(String[] args) throws Exception
	{
		if ( args.length != 2 ) 
		{
			System.out.println( "Usage ConservationSum inAlignment outFile" );
			return;
		}	
		
		Alignment a = new Alignment("1", new File( args[0]), false );
		
		ConservationSum cSum = new ConservationSum(a);
		
		BufferedWriter writer = new BufferedWriter( new FileWriter( new File(
						args[1] )));
		
		writer.write( "i\tj\tscore\n");
		
		for ( int i =0; i < a.getNumColumnsInAlignment(); i++ ) 
			if ( a.columnHasValidResidue(i) ) 
				for ( int j = i + 1; j < a.getNumColumnsInAlignment(); j++ )
					if ( a.columnHasValidResidue(j) ) 
						writer.write( i + "\t" + j + "\t" + cSum.getScore(a, i, j) + "\n" );				
		
		writer.flush();  writer.close();			
	}
	
	public double getScore( int i ) throws Exception
	{
		return absoluteSequenceEntropy[i];
	}
	
	public ConservationSum(Alignment a, String id) throws Exception
	{
		this.id = id;
		this.alignment = a;
		int[][] counts  = a.getCounts();
		absoluteSequenceEntropy = new double[ a.getNumColumnsInAlignment()];
		
		for ( int x=0; x < a.getNumColumnsInAlignment(); x++ )
		{
			absoluteSequenceEntropy[x] = 0.0;
			
			long total = 0;
			
			for ( int y=0; y < MapResiduesToIndex.NUM_VALID_RESIDUES; y++ ) 
				total += counts[x][y];
			
			for ( int y=0; y< MapResiduesToIndex.NUM_VALID_RESIDUES; y++ ) 
			{
				float frequency = ((float) counts[x][y]) / total;
				
				if ( frequency > 0.0 ) 
					absoluteSequenceEntropy[x] += ( frequency * Math.log( frequency)  );
			}
			
			absoluteSequenceEntropy[x] = - absoluteSequenceEntropy[x];
		}

	}
	
	public String getAnalysisName()
	{
		return id;
	}
	
	public boolean isSymmetrical()
	{
		return true;
	}
	
	public boolean reverseSort()
	{
		return true;
	}
	
	public ConservationSum(Alignment a) throws Exception
	{
		this( a, "ConservationSum");
	}
}
