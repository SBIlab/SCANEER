
package covariance.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.StringTokenizer;

import covariance.datacontainers.Alignment;
import covariance.datacontainers.PdbChain;
import covariance.datacontainers.PdbFileWrapper;
import covariance.datacontainers.PdbResidue;
import covariance.parsers.ClustalAlignment;
import covariance.parsers.PFamPdbAnnotationParser;

public class SequenceUtils
{
		private static float[] maxSurfaceAreas = { 
		188.789f,
		200.211f,
		226.145f,
		265.54f,
		301.605f,
		155.998f,
		246.736f,
		247.242f,
		291.399f,
		254.034f,
		280.255f,
		238.508f,
		216.448f,
		252.54f,
		309.349f,
		196.96f,
		215.754f,
		223.66f,
		281.692f,
		300.091f };
					
	// static access only
	private SequenceUtils()
	{
	}
	
	public static void main(String[] args) throws Exception
	{
		for ( int x=0; x< MapResiduesToIndex.NUM_VALID_RESIDUES; x++ ) 
		{
			System.out.println( MapResiduesToIndex.getChar(x) + " " + getMaxSurfaceArea(x) );
		}
	}
	
	public static float getMaxSurfaceArea(int index)
	{
		return maxSurfaceAreas[index];	
	}
	
	public static boolean isHydrophicChar( char inChar ) 
	{
		if ( inChar == 'A' || inChar == 'V' || inChar == 'L' 
				|| inChar == 'I' || inChar =='F' || inChar =='W' 
				|| inChar == 'M' || inChar == 'P') 
					return true;
					
		return false;
	}
	
	public static String stripSpaces( String inString ) 
	{
		StringBuffer buff = new StringBuffer();
		
		StringTokenizer sToken = new StringTokenizer( inString );
		
		while ( sToken.hasMoreElements() )
			buff.append(sToken.nextToken());
		
		return buff.toString();	
	}
	
	public static boolean pairIsHydrophobic( char inChar1, char inChar2 ) 
	{
		if ( isHydrophicChar(inChar1) && isHydrophicChar(inChar2) ) 
			return true;
			
		return false;	
	}
	
	public static String buildFasta(String sequence1, String sequence2 ) 
	{
		StringBuffer buff = new StringBuffer();
		
		buff.append(">1\n");
		buff.append( sequence1 + "\n" );
		buff.append(">2\n");
		buff.append( sequence2 + "\n" );
		
		
		return buff.toString();	
	}
	
	public static PFamPdbAnnotationParser getBestMatchToPdbs(Alignment a, String sequence) throws Exception
	{
		float highestPercentage = -10000;
		PFamPdbAnnotationParser bestHit = null;
		
		for ( int x=0; x< a.getAnnotationParsers().length; x++ )
		{
		
			String pdbPath = ConfigReader.getLocalPdbDirectory() + File.separator + 
								a.getAnnotationParsers()[x].getFourCharId();
			
			PdbFileWrapper pdbFileWrapper = null;
			
			try
			{
				System.out.println("Trying " + a.getAnnotationParsers()[x].getFourCharId());
				pdbFileWrapper = new PdbFileWrapper(a.getAnnotationParsers()[x].getFourCharId());
			} 
			catch(Exception e)
			{
				System.out.println("Could not parse " + a.getAnnotationParsers()[x].getFourCharId() );
			}
			
			if ( pdbFileWrapper != null )
			{
			
				String pdbFragment = SequenceUtils.getPdbStringFragment( a.getAnnotationParsers()[x], pdbFileWrapper );
				ClustalAlignment cAlingment = 
						ClustalWrapper.getClustalAlignment(SequenceUtils.buildFasta(sequence, pdbFragment));
				
				float pairwiseIdentity = cAlingment.getPairwiseIdentity(0,1);
									
				if (  pairwiseIdentity > highestPercentage ) 
				{
					highestPercentage = pairwiseIdentity ;
					bestHit = a.getAnnotationParsers()[x];
					System.out.println("Setting highest to " + pairwiseIdentity );
				}
				
			}
		}
		
		return bestHit;
	}
	
	public static String fileToString( File file ) throws Exception
	{
		BufferedReader reader = new BufferedReader( new FileReader( file ));
		StringWriter writer = new StringWriter();
		
		char[] c = new char[2048];
        int bytesRead;
        
        while ( ( bytesRead = reader.read(c,0,c.length)) != -1 )
			writer.write( c, 0, bytesRead );
		
		return writer.toString();
	}
	
	public static int threeToIndex( String inString ) throws Exception
	{
		return MapResiduesToIndex.getIndex( threeToOne( inString ) );
	}
	
	
	public static String threeToOne( String inString ) throws Exception
	{
		inString = inString.toUpperCase();
		
		if ( inString.equals("GLY") ) 
			return "G";
		
		if ( inString.equals("ALA" )) 
			return "A";
		
		if ( inString.equals("VAL")) 
			return "V";
		
		if ( inString.equals("LEU")) 
			return "L";
		
		if ( inString.equals("ILE" )) 
			return "I";
		
		if ( inString.equals("PHE" )) 
			return "F";
		
		if ( inString.equals("TYR" )) 
			return "Y";
		
		if ( inString.equals("TRP" )) 
			 return "W";
		
		if ( inString.equals("SER" )) 
			return "S";
		
		if ( inString.equals("THR" )) 
			return "T";
		
		if ( inString.equals("MET" )) 
			return "M";
		
		if ( inString.equals("CYS" )) 
			return "C";
		
		if ( inString.equals("ASN" )) 
			return "N";
		
		if ( inString.equals("GLN" )) 
			return "Q";
	
		if ( inString.equals("PRO" ))
			return "P";
		
		if ( inString.equals("ASP" )) 
			return "D";
		
		if ( inString.equals("GLU" ))
			 return "E";
		
		if ( inString.equals("LYS" )) 
			return "K";
		
		if ( inString.equals("ARG" ))
			return "R";
		
		if ( inString.equals("HIS" ))
			return "H";
		
		throw new Exception("Unexpected token '" + inString + "'" );
	}

	public static String getPdbStringFragment( PFamPdbAnnotationParser parser, PdbFileWrapper pdbWrapper )
	{
		StringBuffer buff = new StringBuffer();
		PdbChain chain = (PdbChain) pdbWrapper.getChain(parser.getChainChar());
		
		if ( chain == null ) 
			return null;
		
		for ( int x=parser.getStartPos(); x <= parser.getEndPos(); x++ ) 
		{
			PdbResidue residue = chain.getPdbResidueByPdbPosition(x);
			
			if ( residue != null ) 
				buff.append( residue.getPdbChar());	
		}
		
		return buff.toString();	
	}
	
	public static String getPdbStringFragment( char chainChar, 
												int startPosition, 
												int endPosition, 
												PdbFileWrapper pdbWrapper )
	{
		StringBuffer buff = new StringBuffer();
		PdbChain chain = (PdbChain) pdbWrapper.getChain(chainChar);
		
		if ( chain == null ) 
			return null;
		
		for ( int x=startPosition; x <= endPosition; x++ ) 
		{
			PdbResidue residue = chain.getPdbResidueByPdbPosition(x);
			
			if ( residue != null ) 
				buff.append( residue.getPdbChar());	
		}
		
		return buff.toString();	
	}

	/**  Removes dupliace chains.  If there are two identical chains from the same pdb,
	 *   the shortest one is removed.
	 * 
	 *   Made public and static for testing purposes ( I find I am missing the C++ friend! )
	 */
		public static void trimPdbIDs(List pdbIds) throws Exception
	{
		List listCopy = new ArrayList( pdbIds );
		
		for ( Iterator i = listCopy.iterator();
			   i.hasNext(); ) 
		{
			PFamPdbAnnotationParser s1 = new PFamPdbAnnotationParser( i.next().toString());
			
			for ( Iterator i2 = listCopy.iterator();
				   i2.hasNext(); ) 
			{
				PFamPdbAnnotationParser s2 = new PFamPdbAnnotationParser( i2.next().toString());
				
				if ( ! s1.equals(s2) ) 
				{
					if ( s1.getChainChar()== s2.getChainChar() && 
					     s1.getFourCharId().equals( s2.getFourCharId()) ) 
					{
						if ( s1.getLength() > s2.getLength() ) 
							pdbIds.remove(s2.getPdbAnnotationString());
						else
							pdbIds.remove(s1.getPdbAnnotationString());	
					}
				}	
			}
		}
	}

	

	

}
