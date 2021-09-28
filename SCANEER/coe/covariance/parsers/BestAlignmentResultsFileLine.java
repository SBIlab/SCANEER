package covariance.parsers;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.StringTokenizer;


import covariance.datacontainers.Alignment;
import covariance.datacontainers.AlignmentLine;

public class BestAlignmentResultsFileLine implements Comparable
{
	private String pFamId;
	private String pdbID;
	private char pdbChain;
	private int pdbStartResidue;
	private int pdbEndResidue;
	private float percentIdentity;
	private String alignmentLineId;
	private int numSequencesInAlignment;
	private int alignmentNum;
	
	public static BestAlignmentResultsFileLine getOneFromAlignmentNum(File fileToParse,int alignmentNum) throws Exception
	{
		List list = getBestAlignmentResultsList(fileToParse);
		
		for ( Iterator i = list.iterator();
				i.hasNext(); ) 
		{
			BestAlignmentResultsFileLine rFileLine = (BestAlignmentResultsFileLine) i.next();
			if ( rFileLine.getAlignmentNum() == alignmentNum ) 
				return rFileLine;	
		}
		
		return null;	
	}
	
	
	public static BestAlignmentResultsFileLine getOne( File fileToParse, String pFamID ) throws Exception
	{
		List list = getBestAlignmentResultsList(fileToParse);
		
		for ( Iterator i = list.iterator();
				i.hasNext(); ) 
		{
			BestAlignmentResultsFileLine rFileLine = (BestAlignmentResultsFileLine) i.next();
			if ( rFileLine.getPFamId().equals(pFamID) ) 
				return rFileLine;	
		}
		
		return null;	
	}
	
	public BestAlignmentResultsFileLine( Alignment a,
										  AlignmentLine theLine,
										  PFamPdbAnnotationParser annotationParser,
										  float percentIdentity ) 
	{
		this.pFamId = a.getAligmentID();
		this.pdbID = annotationParser.getFourCharId();
		this.pdbChain = annotationParser.getChainChar();
		this.pdbStartResidue = annotationParser.getStartPos();
		this.pdbEndResidue = annotationParser.getEndPos();
		this.percentIdentity = percentIdentity;	
		this.alignmentLineId = theLine.getIdentifier();
		this.numSequencesInAlignment = a.getNumSequencesInAlignment();
	}	
	
	public BestAlignmentResultsFileLine( Alignment a,
										  AlignmentLine theLine,
										  String fourCharPDBId,
										  char chainChar,
										  int pdbStart,
										  int pdbEnd ) 
	{
		this.pFamId = a.getAligmentID();
		this.pdbID = fourCharPDBId;
		this.pdbChain = chainChar;
		this.pdbStartResidue = pdbStart;
		this.pdbEndResidue = pdbEnd;
		this.alignmentLineId = theLine.getIdentifier();
		this.numSequencesInAlignment = a.getNumSequencesInAlignment();
			
	}											  							  
	
	public BestAlignmentResultsFileLine( String fileLine ) throws Exception
	{
		StringTokenizer sToken = new StringTokenizer( fileLine, "\t" );
		this.pFamId= sToken.nextToken();
		this.alignmentNum = Integer.parseInt( sToken.nextToken()); 
		sToken.nextToken(); // numColumnsInAlignment
		this.numSequencesInAlignment = Integer.parseInt( sToken.nextToken()); // numSequencesInAlignment
		
		if ( sToken.hasMoreElements() ) 
		{
			this.alignmentLineId = sToken.nextToken();
			sToken.nextToken(); // ungapped Residue Length
			this.pdbID = sToken.nextToken();
			
			String chainToken = sToken.nextToken();
				
			if ( chainToken == null || chainToken.trim().length() == 0 )
				this.pdbChain = ' ';
			else
				this.pdbChain = chainToken.charAt(0);
			
			this.pdbStartResidue = Integer.parseInt( sToken.nextToken());  // pdbStartResidue
			this.pdbEndResidue = Integer.parseInt( sToken.nextToken()); //pdb End Reisdue
			
			sToken.nextToken();  // pdbLength
			
			this.percentIdentity = Float.parseFloat( sToken.nextToken() );
			sToken.nextToken();  // minimun length
			
			if ( sToken.hasMoreTokens() ) 
				throw new Exception("Parsing error");		
		}
		
	}

	public String getAlignmentLineId()
	{
		return alignmentLineId;
	}

	public char getPdbChain()
	{
		return pdbChain;
	}

	public int getPdbLength()
	{
		return pdbEndResidue - pdbStartResidue;
	}

	public int getPdbEndResidue()
	{
		return pdbEndResidue;
	}

	public String getPdbID()
	{
		return pdbID;
	}

	public int getPdbStartResidue()
	{
		return pdbStartResidue;
	}

	public float getPercentIdentity()
	{
		return percentIdentity;
	}

	public String getPFamId()
	{
		return pFamId;
	}
	
	public static List<BestAlignmentResultsFileLine> getBestAlignmentResultsList(String filePath)
		throws Exception
	{
		return getBestAlignmentResultsList(new File(filePath));
	}
	
	public static HashSet<String> getPfamIdsAsSet(List<BestAlignmentResultsFileLine> list) throws Exception
	{
		HashSet<String> set = new HashSet<String>();
		
		for ( BestAlignmentResultsFileLine rFileLine : list )
			set.add(rFileLine.getPFamId());
		
		return set;
	}
	
	public static List<BestAlignmentResultsFileLine> getBestAlignmentResultsList(File fileToParse) throws Exception
	{
		List<BestAlignmentResultsFileLine> list = new ArrayList<BestAlignmentResultsFileLine>();
		
		System.out.println("Opening "  + fileToParse.getAbsolutePath() );
		BufferedReader reader = new BufferedReader( new FileReader( fileToParse));
		
		reader.readLine();
		
		String nextLine = reader.readLine();
		
		while ( nextLine != null && nextLine.trim().length() != 0 ) 
		{
			list.add( new BestAlignmentResultsFileLine( nextLine ));
			nextLine = reader.readLine();
		}	
		
		return list;
	}
	
	public static BestAlignmentResultsFileLine getOneById(List<BestAlignmentResultsFileLine> list,
			String id) throws Exception
	{
		for ( BestAlignmentResultsFileLine rFileLine : list )
			if ( rFileLine.getPFamId().equals(id))
				return rFileLine;
			
		throw new Exception("Could not find" + id);
	}
	
	public static BestAlignmentResultsFileLine getOneByID( File fileToParse, String pFamID ) throws Exception
	{
		for ( Iterator i = getBestAlignmentResultsList(fileToParse).iterator();
					i.hasNext(); ) 
		{
			BestAlignmentResultsFileLine rFileLine = (BestAlignmentResultsFileLine) i.next();
			
			if ( rFileLine.getAlignmentLineId() != null && rFileLine.getPFamId().equals(pFamID) ) 
				return rFileLine;
		}
		
		return null;	
	}
	
	public int getNumSequencesInAlignment() 
	{
		return numSequencesInAlignment;
	}

	public int getAlignmentNum()
	{
		return alignmentNum;
	}

	public int compareTo(Object o)
	{
		return this.pFamId.toUpperCase().compareTo( ((BestAlignmentResultsFileLine)o).pFamId.toUpperCase() );
	}

	public void setPdbID(String pdbID)
	{
		this.pdbID = pdbID;
	}

	public void setPdbChain(char pdbChain)
	{
		this.pdbChain = pdbChain;
	}

}
