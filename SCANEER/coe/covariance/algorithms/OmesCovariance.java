package covariance.algorithms;

import java.io.*;
import java.util.*;

import covariance.datacontainers.*;
import covariance.utils.*;

public class OmesCovariance implements ScoreGenerator
{ 
	private String idString;
	private boolean tossSuspiciousPoints = false;
	private Alignment a;
	
	public OmesCovariance( Alignment a, String idString ) 
	{
		this.idString = idString;
		this.a = a;
	}
	
	public static void main(String[] args) throws Exception
	{
		if ( args.length != 2 ) 
		{
			System.out.println( "Usage OmesCovariance inAlignment outFile" );
			return;
		}	
		
		Alignment a = new Alignment("1", new File( args[0]), false );
		
		OmesCovariance omes = new OmesCovariance(a);
		
		BufferedWriter writer = new BufferedWriter( new FileWriter( new File(
						args[1] )));
		
		writer.write( "i\tj\tscore\n");
		
		for ( int i =0; i < a.getNumColumnsInAlignment(); i++ ) 
			if ( a.columnHasValidResidue(i) ) 
				for ( int j = i + 1; j < a.getNumColumnsInAlignment(); j++ )
					if ( a.columnHasValidResidue(j) ) 
					writer.write( i + "\t" + j + "\t" + omes.getScore(a, i, j) + "\n" );				
		
		writer.flush();  writer.close();			
	}
	
	public OmesCovariance(Alignment a)
	{
		this(a, "OMES");
	}
	
	private static class PairsClass
	{
		String pair;
		int num;
		
		PairsClass( String pair ) 
		{
			this.pair = pair;
			num = 0;
		}
		
		void increment()
		{
			num++;
		}
	}
	
	private static Random random = new Random();
	
	public String getAnalysisName()
	{
		return idString;
	}
	
	/**  The key is the two char string representing that pair
	 *   The values is a PairsClass object
	 */
	public HashMap getPairs( String iString, 
							 String jString,
							 int[] iFreqs,
							 int[] jFreqs) throws Exception
	{
		HashMap pairs = new HashMap();
		
		for ( int i=0; i< iFreqs.length; i++ ) 
			for ( int j=0; j < jFreqs.length; j++ )
				if ( iFreqs[i] > 0 && jFreqs[j] > 0 ) 
				{
					String pair = "" + MapResiduesToIndex.getChar( i ) + MapResiduesToIndex.getChar( j );
					pairs.put( pair, new PairsClass( pair ));
				}											
		
		for ( int x=0; x< iString.length(); x++ ) 
		{
			char iChar = iString.charAt(x);
			char jChar = jString.charAt(x);
			
			String pair = "" + iChar+ jChar;
			((PairsClass) pairs.get( pair )).increment();
		}
		
		return pairs;
	}
	
	public static int[] getFrequencies(String string) throws Exception
	{
		int[] frequencies = new int[ MapResiduesToIndex.NUM_VALID_RESIDUES ];
		
		for ( int x=0; x< string.length(); x++ ) 
		{
			frequencies[ MapResiduesToIndex.getIndex( string.charAt(x) )]++;
		}
		
		return frequencies;
	}
	
	public double getCovarianceScore( HashMap pairs,
									  int numValidSequences,
									  int[] iFrequencies,
									  int[] jFrequencies,
									  int i,
									  int j) throws Exception
	{
		if ( numValidSequences <= 0 )
			return 0;
		double covarianceScore = 0;
		
		for ( Iterator it = pairs.values().iterator();
			  it.hasNext(); ) 
		{
			PairsClass pClass= (PairsClass) it.next();
			char iChar = pClass.pair.charAt(0);
			char jChar = pClass.pair.charAt(1);
			
			// the frequency of residue x at position i
			float fxi = ((float)iFrequencies[ MapResiduesToIndex.getIndex(iChar)]) / numValidSequences;
			
			//System.out.println( "Freq of " + MapResiduesToIndex.getChar( cvh.getIResidue()) + "=" + fxi );
			
			// the frequency of residue y at poistion j
			float fyj = ((float) jFrequencies[ MapResiduesToIndex.getIndex(jChar)] ) / numValidSequences;
			
			// expected number of sequences that cotain an x at pos i and a y at position j
			float nex = numValidSequences * fxi * fyj;

			float score = pClass.num - nex;
			
			covarianceScore += score * score;
			//System.out.println( iChar + "\t" + jChar + "\t" + pClass.num + "\t" + nex + "\t" + score*score );
		}
		
		return covarianceScore / numValidSequences;
	}
	
	
	public List getSubScores( Alignment alignment, int i, int j ) throws Exception
	{
		if ( alignment != this.a )
			throw new Exception("Please call on the same alignment used in the constructor");
		
		String iString = alignment.getColumnAsString(i);
		String jString= alignment.getColumnAsString(j);
		
		if ( iString.length() != jString.length() )
			throw new Exception("Logic error");
		
		StringBuffer iStringNew = new StringBuffer();
		StringBuffer jStringNew= new StringBuffer();
		
		for ( int x=0; x< iString.length(); x++ ) 
		{
			char iListChar = iString.charAt(x);;
			char jListChar = jString.charAt(x);
			
			if ( MapResiduesToIndex.isValidResidueChar( iListChar ) &&
				 MapResiduesToIndex.isValidResidueChar( jListChar ) )
			{
				iStringNew.append( iListChar );
				jStringNew.append( jListChar );
			}
		}
		
		iString = iStringNew.toString();
		jString = jStringNew.toString();
		int[] iFrequencies = getFrequencies( iString );
		int[] jFrequencies = getFrequencies( jString);
		
		HashMap pairs = getPairs( iString, jString, iFrequencies, jFrequencies );
		int numValidSequences = iString.length();
		
		List list = new ArrayList();
		
		for ( Iterator it = pairs.values().iterator();
			  it.hasNext(); ) 
		{
			PairsClass pClass= (PairsClass) it.next();
			char iChar = pClass.pair.charAt(0);
			char jChar = pClass.pair.charAt(1);
			
			// the frequency of residue x at position i
			float fxi = ((float)iFrequencies[ MapResiduesToIndex.getIndex(iChar)]) / numValidSequences;
			
			//System.out.println( "Freq of " + MapResiduesToIndex.getChar( cvh.getIResidue()) + "=" + fxi );
			
			// the frequency of residue y at poistion j
			float fyj = ((float) jFrequencies[ MapResiduesToIndex.getIndex(jChar)] ) / numValidSequences;
			
			// expected number of sequences that cotain an x at pos i and a y at position j
			float nex = numValidSequences * fxi * fyj;

			float score = pClass.num - nex;
			
			list.add( new AlignmentSubScore(iChar, jChar, pClass.num, nex, score * score ));
		}
		
		return list;
	}
	
	public double getScore( Alignment alignment, int i, int j ) throws Exception
	{
		if ( alignment != this.a )
			throw new Exception("Please call on the same alignment used in the constructor");
		
		String iString = alignment.getColumnAsString(i);
		String jString= alignment.getColumnAsString(j);
		
		if ( iString.length() != jString.length() )
			throw new Exception("Logic error");
		
		StringBuffer iStringNew = new StringBuffer();
		StringBuffer jStringNew= new StringBuffer();
		
		for ( int x=0; x< iString.length(); x++ ) 
		{
			char iListChar = iString.charAt(x);;
			char jListChar = jString.charAt(x);
			
			if ( MapResiduesToIndex.isValidResidueChar( iListChar ) &&
				 MapResiduesToIndex.isValidResidueChar( jListChar ) )
			{
				iStringNew.append( iListChar );
				jStringNew.append( jListChar );
			}
		}
		
		iString = iStringNew.toString();
		jString = jStringNew.toString();
		int[] iFrequencies = getFrequencies( iString );
		int[] jFrequencies = getFrequencies( jString);
		
		HashMap pairs = getPairs( iString, jString, iFrequencies, jFrequencies );
		return getCovarianceScore( pairs, iString.length(), iFrequencies, jFrequencies, i, j );
		
	}
	
	public boolean isSymmetrical()
	{
		return true;
	}
	
	public boolean reverseSort()
	{
		return false;
	}
	public void setTossSuspiciousPoints(boolean tossSuspiciousPoints)
	{
		this.tossSuspiciousPoints = tossSuspiciousPoints;
	}

	public boolean isTossSuspiciousPoints()
	{
		return tossSuspiciousPoints;
	}
}
