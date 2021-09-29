package covariance.algorithms;

import java.io.*;
import java.util.*;

import covariance.datacontainers.*;

public class RandomScore implements ScoreGenerator
{
	Random random = new Random();
	String idString;
	
	public static void main(String[] args) throws Exception
	{
		if ( args.length != 2 ) 
		{
			System.out.println( "Usage RandomScore inAlignment outFile" );
			return;
		}	
		
		Alignment a = new Alignment("1", new File( args[0]), false );
		
		RandomScore rScore = new RandomScore();
		
		BufferedWriter writer = new BufferedWriter( new FileWriter( new File(
						args[1] )));
		
		writer.write( "i\tj\tscore\n");
		
		for ( int i =0; i < a.getNumColumnsInAlignment(); i++ ) 
			if ( a.columnHasValidResidue(i) ) 
				for ( int j = i + 1; j < a.getNumColumnsInAlignment(); j++ )
					if ( a.columnHasValidResidue(j) ) 
						writer.write( i + "\t" + j + "\t" + rScore.getScore(a, i, j) + "\n" );				
		
		writer.flush();  writer.close();			
	}
	
	public double getScore( Alignment a, int i, int j ) throws Exception
	{
		return random.nextDouble();
	}
	
	public RandomScore()
	{
		this("random");
	}
	
	public RandomScore(String idString)
	{
		this.idString = idString;
	}
	
	public String getAnalysisName()
	{
		return idString;
	}
	
	public boolean isSymmetrical()
	{
		return true;
	}
	
	public boolean reverseSort()
	{
		return true;
	}
}
