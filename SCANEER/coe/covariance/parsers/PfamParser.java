package covariance.parsers;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

import covariance.datacontainers.Alignment;
import covariance.datacontainers.AlignmentLine;
import covariance.utils.ConfigReader;
import covariance.utils.SequenceUtils;

public class PfamParser
{
	private String lastLineRead = null;
	private boolean finished = false;
	private BufferedReader pFamReader;
	private int numLinesRead=0;
	
	public boolean isFinished()
	{
		return finished;
	}
	
	public PfamParser() throws Exception
	{
		pFamReader = new BufferedReader( new FileReader( new File( ConfigReader.getFullPfamPath())));	
	}
	
	private AlignmentLine getAlignmentLine() throws Exception
	{
		//System.out.println("Got initial line of " + lastLineRead );
		
		int slashIndex = lastLineRead.indexOf("/");
		
		if ( lastLineRead.indexOf("/") == -1 ) 
			throw new Exception("Error!  Expecting a forward slash at line" + numLinesRead );
		
		String id = lastLineRead.substring(0, slashIndex );
		
		String firstLine = lastLineRead.substring( slashIndex + 1 );
		int lowerBound = Integer.parseInt( firstLine.substring(0, firstLine.indexOf("-") ));
		firstLine =  firstLine.substring( firstLine.indexOf("-"));
		StringTokenizer sToken = new StringTokenizer( firstLine );
		int upperBound = Integer.parseInt( sToken.nextToken() );
		
		AlignmentLine aLine = new AlignmentLine( id, sToken.nextToken() );
		
		return aLine;	
	}


	public static Alignment getOneAlignment( String alignmentName ) throws Exception
	{
		PfamParser pFamParser = new PfamParser();
		
		while( ! pFamParser.isFinished() )
		{
			Alignment a = pFamParser.getNextAlignment();
			System.out.println("Scanned " + a.getAligmentID() );
			
			if ( a.getAligmentID().equals( alignmentName) )
				return a;	
		}
		
		return null;
	}

	/**  gets the next alignment from the reader.
	 * 
	 *   The cursor in the reader should be set to a line that starts with "# STOCKHOLM 1.0"
	 *   or this will throw.
	 * 
	 *   Returns null when there are no more alignments to read ( after which isFinished() will return true)
	 */
	public Alignment getNextAlignment() throws Exception
	{
		List pdbIds = new ArrayList();
		lastLineRead = getNextLine();
		
		if ( lastLineRead == null ) 
			return null;
		
		if ( ! lastLineRead.trim().equals("# STOCKHOLM 1.0") )
			throw new Exception("Error!  First line should be # STOCKHOLM 1.0 at "+ numLinesRead );
		
		lastLineRead = getNextLine();
		
		if ( ! lastLineRead.startsWith("#=GF ID" )) 
			throw new Exception("Error!  Second line should start #=GF ID at "+ numLinesRead );
		
		int index = lastLineRead.indexOf("#=GF ID" );
		
		if ( index == -1 ) 
			throw new Exception("Logic error");
		
		String id = lastLineRead.substring(7).trim();
		
		while ( lastLineRead.startsWith("#")) 
		{
			StringTokenizer sToken = new StringTokenizer(lastLineRead);
			if ( sToken.nextToken().equals("#=GF") &&
				 sToken.nextToken().equals("DR") &&
				 sToken.nextToken().startsWith("PDB")  ) 
				pdbIds.add( lastLineRead );
			
			lastLineRead = getNextLine();
		}
		
		List alignmentLines = new ArrayList();
		while ( lastLineRead != null && ! lastLineRead.trim().equals("//") ) 
		{
			if ( ! lastLineRead.startsWith("#" )) 
				alignmentLines.add( getAlignmentLine() );
			
			lastLineRead = getNextLine();
		}
		
		Alignment a = new Alignment( id, alignmentLines );
		SequenceUtils.trimPdbIDs(pdbIds);
		a.setPdbIds( pdbIds);
		
		return a;
	}
	
	private String getNextLine() throws Exception
	{
		if ( finished ) 
			throw new Exception("Error!  Nothing left to read");
		
		numLinesRead++;
		String nextLine = pFamReader.readLine();
		if ( nextLine == null || nextLine.trim().length() == 0 )
			finished = true;
			
		//System.out.println("Reading " + nextLine );
		return nextLine;
	}
}
