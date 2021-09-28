package covariance.algorithms;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.StringTokenizer;

import covariance.datacontainers.Alignment;
import covariance.datacontainers.AlignmentLine;
import covariance.utils.ConfigReader;
import covariance.utils.MapResiduesToIndex;

public class McBASCCovariance implements ScoreGenerator
{	
	private static int[][] maxhomMetric;
	
	private float[] averages;
	private float[] sds;
	private Alignment a;
	private int numSquared;
	
	public McBASCCovariance( Alignment a) throws Exception
	{
		getMaxhomMetric();  // initalize array
		this.a = a;	
		this.numSquared = a.getNumSequencesInAlignment() * a.getNumSequencesInAlignment();
		
		averages = new float[a.getNumColumnsInAlignment()];
		sds = new float[a.getNumColumnsInAlignment()];
		
		for ( int i=0; i < a.getNumColumnsInAlignment(); i++ ) 
		{
			long sum = 0;
			long sumSquared = 0;
			int numComparisons = 0;
			
			for ( int x=0; x < a.getNumSequencesInAlignment(); x++ ) 
			{
				String lineX = ((AlignmentLine) a.getAlignmentLines().get(x)).getSequence();
			
				char xi = lineX.charAt(i);
				
				// could make the algorithm 2X faster by having this loop be 
				// "for ( int y=x+1; y < a.getNumSequencesInAlignment(); y++ ) "
				for ( int y=0; y < a.getNumSequencesInAlignment(); y++ ) 
				{
					
					// should probably uncomment this so columns are not compared to themselves
					// but, for better or worse, that's not the way we described the algorithm in the paper...
					//if ( y != x ) 
					{
						String lineY = ((AlignmentLine) a.getAlignmentLines().get(y)).getSequence();;
						
						char yi = lineY.charAt(i);
						
						if ( MapResiduesToIndex.isValidResidueChar(xi) &&
						MapResiduesToIndex.isValidResidueChar(yi) ) 
						{
							numComparisons++;
							int metric = 
								maxhomMetric[ MapResiduesToIndex.getIndex(xi) ][MapResiduesToIndex.getIndex(yi)];
							sum+= metric;
							sumSquared += metric * metric;	
						}	
					}
				}
				
			}
			
			if ( numComparisons == 0 ) 
			{
				averages[i] =0;
				sds[i] = 0;
			}
			else if ( numComparisons == 1 ) 
			{
				averages[i] = sum;
				sds[i] = 0;
			}
			else
			{	
				averages[i]= ((float) sum) / numComparisons;
				float top = sumSquared;
				top -= ((float)sum) * sum / numComparisons;
				top /= ((float)numComparisons) - 1;
				
				if ( Math.abs(top) < 0.0000000001 ) 
					sds[i] = 0;
				else
					sds[i] = (float) Math.sqrt(top);	
					
				if ( Float.isNaN(sds[i])|| Float.isInfinite(sds[i]) )
					sds[i] = 0;
			}
		}				
	}
	
	public synchronized static int[][] getMaxhomMetric() throws Exception
	{
		if ( maxhomMetric == null ) 
		{
			maxhomMetric = new int[MapResiduesToIndex.NUM_VALID_RESIDUES][MapResiduesToIndex.NUM_VALID_RESIDUES];
		
			for ( int x=0; x< MapResiduesToIndex.NUM_VALID_RESIDUES; x++ ) 
				for ( int y=0; y < MapResiduesToIndex.NUM_VALID_RESIDUES; y++) 
					maxhomMetric[x][y] = -99;
					
			fillMatrixFromFile(maxhomMetric);
			
			for ( int x=0; x< MapResiduesToIndex.NUM_VALID_RESIDUES; x++ ) 
				for ( int y=0; y < MapResiduesToIndex.NUM_VALID_RESIDUES; y++) 
					if ( maxhomMetric[x][y] == -99 ) 
						throw new Exception("Unfilled positon " + x + " " + y );
						
		}
		
		return maxhomMetric;
	}
	
	private static void fillMatrixFromFile(int[][] matrix) throws Exception
	{
		BufferedReader reader = new BufferedReader( new FileReader( new File(
				ConfigReader.getHomeDirectory() + File.separator + "data" + File.separator +
						"Maxhom_McLachlan.metric" )));
		
		String firstLine = reader.readLine();
		StringTokenizer sToken = new StringTokenizer( firstLine );
		char[] topRow = new char[MapResiduesToIndex.NUM_VALID_RESIDUES];
						
		for ( int x=0; x < MapResiduesToIndex.NUM_VALID_RESIDUES; x++ ) 
		{
			String nextChar = sToken.nextToken();
			
			if ( nextChar.length() != 1 ) 
				throw new Exception("Unexpected token " + nextChar );
				
			topRow[x] = nextChar.charAt(0);
		}
		
		for ( int x=0; x < MapResiduesToIndex.NUM_VALID_RESIDUES; x++ ) 
		{
			StringTokenizer lineTokens = new StringTokenizer( reader.readLine() );
			
			String nextChar = lineTokens.nextToken();
			
			if ( nextChar.length() != 1 ) 
				throw new Exception("Unexpected token " + nextChar );
					
			for ( int y=0; y < MapResiduesToIndex.NUM_VALID_RESIDUES; y++ ) 
			{
				matrix[ MapResiduesToIndex.getIndex( nextChar.charAt(0))][ MapResiduesToIndex.getIndex( topRow[y])] = 
							(int) Float.parseFloat(lineTokens.nextToken());
			}
		}
	}
	
	public String getAnalysisName()
	{
		return "McBASC";
	}
	
	public static void main(String[] args) throws Exception
	{
		if ( args.length != 2 ) 
		{
			System.out.println( "Usage McBASCCovariance inAlignment outFile" );
			return;
		}	
		
		Alignment a = new Alignment("1", new File( args[0]), false );
		a = a.getFilteredAlignment(90);
		
		McBASCCovariance mcBasc  = new McBASCCovariance(a);
		
		BufferedWriter writer = new BufferedWriter( new FileWriter( new File(
						args[1] )));
		
		writer.write( "i\tj\tscore\n");
		
		for ( int i =0; i < a.getNumColumnsInAlignment(); i++ ) 
				if ( a.columnHasValidResidue(i) ) 
					for ( int j = i + 1; j < a.getNumColumnsInAlignment(); j++ )
						if ( a.columnHasValidResidue(j) ) 
							writer.write( i + "\t" + j + "\t" + mcBasc.getScore(a, i, j) + "\n" );				
			
		writer.flush();  writer.close();			
	}
	
	public float getMatrixTopScore( int i, int j, char xi, char xj, char yi, char yj ) 
		throws Exception
	{
	
		if ( MapResiduesToIndex.isValidResidueChar(xi) &&
						MapResiduesToIndex.isValidResidueChar(yi) && 
						 MapResiduesToIndex.isValidResidueChar(xj) && 
						  MapResiduesToIndex.isValidResidueChar(yj) ) 
		{
			int iScore = maxhomMetric[MapResiduesToIndex.getIndex(xi)][MapResiduesToIndex.getIndex(yi)];
			int jScore = maxhomMetric[MapResiduesToIndex.getIndex(xj)][MapResiduesToIndex.getIndex(yj)];
				
			float top = (iScore- averages[i]);		
			top *= (jScore - averages[j]);
			//System.out.println( i + " " + j + " " + xi + " " + xj + " " + yi + " " + yj + " " + top );
			return top;
		}
	
		return 0;	
	}

	/*
	 This method makes somes approximations for performance which don't really have to be made.
	 These approximations effect alignments with few sequences in them ( see comments below).
	 A cleaner version could return exact values in less time.
	 But McBASC was working pretty well in our paper (in which we only used alignments with >50 sequences
	 so the approximations we used were nearly exact), so I basically didn't want to mess with it.
	 As always, your milage may vary.  If you are using smaller alignments, you can fix these approximations
	 as indicated in the comments in this method.
	 */
	public double getScore(Alignment a, int i, int j) throws Exception
	{
		if ( a != this.a ) 
			throw new Exception("Please call on alignment passed into constructor!");
			
		// if either column is perfectly conserved, move it somewhere out of the way
		if ( sds[i] == 0 || sds[j] == 0 ) 
			return NO_SCORE;
			
		float sum =0;
		
		for ( int x=0; x < a.getNumSequencesInAlignment(); x++ ) 
		{
			String lineX = ((AlignmentLine) a.getAlignmentLines().get(x)).getSequence();
			
			char xi = lineX.charAt(i);
			char xj = lineX.charAt(j);
			
			// a small bug here.  When calculating the means and standard deviations in 
			// the constructor, we allowed each residue to be compared to itself
			// here, when calcuating the final (i,j) score, we don't compare residues to themselves
			// for any alignment with more than a non-trivial number of sequences ( > ~5-10), this bug
			// makes essentially no difference
			for ( int y=x + 1; y < a.getNumSequencesInAlignment(); y++ ) 
			{
				String lineY = ((AlignmentLine) a.getAlignmentLines().get(y)).getSequence();
					
				char yi = lineY.charAt(i);
				char yj = lineY.charAt(j);
				
				sum += getMatrixTopScore(i, j, xi, xj, yi, yj);
			}
		}
		
		// the factor of 2 multiplication is an approximation for performance.
		// the half triangle in the matrix does not have exactly half the comparisons of the
		// full triangle.  This approximation will approach being exact 
		// as the number of sequences in the alignment increases.
		// For small alignments with few sequences, however, this approximation may effect the score
		// Looking at this now, it seems kind of sloppy, but because we only used large aligments in the paper
		// ( >50 sequences), it made almost no difference to our results
		return 2 * Math.abs( sum / ( this.numSquared * sds[i] * sds[j]));
	}

	public boolean isSymmetrical()
	{
		return true;
	}

	public boolean reverseSort()
	{
		return false;
	}
	
	/** getters for testing
	 */
	public float[] getAverages()
	{
		return averages;
	}

	public float[] getSds()
	{
		return sds;
	}
	public int getNumSquared()
	{
		return numSquared;
	}
}
