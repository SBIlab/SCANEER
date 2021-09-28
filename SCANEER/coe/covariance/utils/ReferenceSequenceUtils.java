package covariance.utils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;

import scripts.PdbDownload;
import covariance.datacontainers.Alignment;
import covariance.datacontainers.AlignmentLine;
import covariance.datacontainers.AlignmentPdbMatch;
import covariance.datacontainers.ClustalLine;
import covariance.datacontainers.PdbChain;
import covariance.datacontainers.PdbFileWrapper;
import covariance.datacontainers.PdbResidue;
import covariance.datacontainers.ReferenceSequence;
import covariance.datacontainers.ReferenceSequenceResidue;
import covariance.parsers.BestAlignmentResultsFileLine;
import covariance.parsers.ClustalAlignment;
import covariance.parsers.PFamPdbAnnotationParser;

public class ReferenceSequenceUtils
{
	public static final String ALIGNMENT_LINE_WORD="d_*_*_12";
	public static final char CHAIN_SPACE_CHAR='-';
	public static float IDENTITY_CUTOFF = 95f;
	
	

	public static String buildFasta( Alignment a,
										AlignmentLine theLine,
										BestAlignmentResultsFileLine rFileLine,
										PdbFileWrapper pdbWrapper ) throws Exception
	{
		StringBuffer buff = new StringBuffer();
		
		buff.append( ">" + ALIGNMENT_LINE_WORD+ "\n" );
		buff.append( theLine.getUngappedSequence() + "\n" );
		
		boolean foundAtLeastOne = false;
		
		for ( int x=0; x< a.getAnnotationParsers().length; x++) 
		{
			PFamPdbAnnotationParser parser= a.getAnnotationParsers()[x];
			if ( parser.getFourCharId().equals(rFileLine.getPdbID() ) )
			{
				foundAtLeastOne = true;
				String pdbFrag = SequenceUtils.getPdbStringFragment( parser, pdbWrapper);			
				
				char chainChar = parser.getChainChar();
				if ( chainChar == ' ' ) 
					chainChar = CHAIN_SPACE_CHAR;	
					
				buff.append( ">" + chainChar + "\n");	
				buff.append(  pdbFrag + "\n" );
			}
		}
		
		
		if ( ! foundAtLeastOne ) 
			throw new Exception("Could not find any chains for " + rFileLine.getPdbID() );
		
		return buff.toString();
		
	}
	
	public static AlignmentLine getBestAlignmentLine( Alignment a,
													    String sequence ) throws Exception
	{
		float highestVal = -10000;
		
		AlignmentLine aLine = null;
		
		for ( Iterator i = a.getAlignmentLines().iterator();
						i.hasNext(); )
		{
			AlignmentLine aPossibleLine = (AlignmentLine) i.next();
			ClustalAlignment cAlignment = ClustalWrapper.getClustalAlignment( 
									SequenceUtils.buildFasta(sequence, aPossibleLine.getUngappedSequence()) );
									
			float pairwiseIdentity = cAlignment.getPairwiseIdentity(0,1);
			
			if ( pairwiseIdentity > highestVal )
			{
				System.out.println(" got " + aPossibleLine.getIdentifier()  + " " + pairwiseIdentity );
				aLine = aPossibleLine;
				highestVal = pairwiseIdentity;
			}
		}
		
		return aLine;
	}				
	
	
	private static List getClustalLines( Alignment a, 
									 ClustalAlignment clustalAlignment,
									 BestAlignmentResultsFileLine rFileLine )
		throws Exception
	{
		List list = new ArrayList();
		boolean foundAlignmentSequence = false;
		
		for ( int x=0; x< clustalAlignment.getNumSequences(); x++ )
		{
			int startPosition = 0;
			boolean isAlignmentSequence = false;
			
			if ( ! clustalAlignment.getSequenceNames()[x].equals( ALIGNMENT_LINE_WORD ) ) 
			{
				char chainChar = clustalAlignment.getSequenceNames()[x].charAt(0);
			
				if ( chainChar == CHAIN_SPACE_CHAR )
					chainChar = ' ';
				
				PFamPdbAnnotationParser parser = a.getAPdbIdParser(rFileLine.getPdbID(), chainChar);
			
				if ( parser == null ) 
					throw new Exception("Error!  Can not find a chain for " + 
								rFileLine.getPdbID() + " " + chainChar);
								
				startPosition = parser.getStartPos();
				isAlignmentSequence = false;

			}
			else
			{
				foundAlignmentSequence = true;
				isAlignmentSequence = true;
			}
			
			list.add(new ClustalLine( clustalAlignment.getSequenceNames()[x], 
										x, startPosition, isAlignmentSequence));
		}
		
		if ( ! foundAlignmentSequence ) 
			throw new Exception("Error!  Could not find alignmentSequence in " + a.getAligmentID() );
		
		
		return list;
	}

	private static PdbChain[] getPdbChains(  
									 Alignment a, 
									 ClustalAlignment clustalAlignment,
									 BestAlignmentResultsFileLine rFileLine,
									 PdbFileWrapper pdbWrapper ) throws Exception
	{
		PdbChain[] chains = new PdbChain[ clustalAlignment.getNumSequences() - 1];
		
		if ( chains.length == 0 ) 
			throw new Exception("Could not find any chains in clustalAlignment for " + a.getAligmentID() );
		
		int numDone = 0;
		for ( int x=0; x< clustalAlignment.getNumSequences(); x++ )
		{
			if ( ! clustalAlignment.getSequenceNames()[x].equals(ALIGNMENT_LINE_WORD))
			{
 				char chainChar = clustalAlignment.getSequenceNames()[x].charAt(0);
			
				if ( chainChar == CHAIN_SPACE_CHAR )
					chainChar = ' ';
			
				chains[numDone] = pdbWrapper.getChain(chainChar);
			
				if ( chains[numDone] == null ) 
					throw new Exception("Could not find a chain for " + a.getAligmentID() + " " + chainChar );
					
				numDone++;
			}
		}
		
		return chains;

	}
	public static AlignmentLine getBestMatchingAlignmentLine(Alignment a,
															  PdbFileWrapper pdbWrapper,
															  char chainChar,
															  int startPosition, 
															  int endPosition ) throws Exception
	{
		String pdbFragment = SequenceUtils.getPdbStringFragment(chainChar,startPosition, endPosition, pdbWrapper );
		float bestSoFar = -1;
		AlignmentLine theLine = null;
		
		for ( Iterator i = a.getAlignmentLines().iterator();
							i.hasNext(); ) 
		{
			AlignmentLine aLine = (AlignmentLine) i.next();
			String aLineUngappedSequence = aLine.getUngappedSequence();
			int alignmentLineLength =  aLineUngappedSequence.length();
					
			String fasta = ReferenceSequenceUtils.getFastaString( aLine, aLineUngappedSequence,
										 pdbWrapper.getFourCharId() + "_" + chainChar, pdbFragment);
			ClustalAlignment clustalAlignment = ClustalWrapper.getClustalAlignment( fasta);					
			float pairwiseIdentity = clustalAlignment.getPairwiseIdentity(0,1);
			
			if ( pairwiseIdentity > bestSoFar )
			{
				theLine = aLine;
				System.out.println("\t\tmatched " + pairwiseIdentity + " for " + theLine.getIdentifier());
				bestSoFar = pairwiseIdentity;
			}
		}
		
		return theLine;		
	}

	private static ReferenceSequenceResidue addFirstLineResidue( 
															Alignment a,
															char[][] alignmentChars,
															ClustalLine clustaledLine,
															int x,
															String cappedLine) throws Exception
	{
		ReferenceSequenceResidue rSeqResidue = null;
		
		if ( MapResiduesToIndex.isValidResidueChar( alignmentChars[clustaledLine.getClustalIndex()][x]) ) 
		{
			rSeqResidue = 
				new ReferenceSequenceResidue( alignmentChars[clustaledLine.getClustalIndex()][x], 
													clustaledLine.getPositionInOriginatingSequence());
			
			if ( cappedLine.charAt(clustaledLine.getPositionInOriginatingSequence()) 
								!= alignmentChars[clustaledLine.getClustalIndex()][x] )
			{
				throw new Exception("Error!  Unequal chars for " +
											a.getAligmentID() + " " 
											+ cappedLine.charAt(clustaledLine.getPositionInOriginatingSequence()) 	
											+ " vs " + alignmentChars[clustaledLine.getClustalIndex()][x] );
			}
				
			clustaledLine.incrementPositionInOriginalSequence();
				
			while(  clustaledLine.getPositionInOriginatingSequence() < cappedLine.length() && 
				!MapResiduesToIndex.isValidResidueChar( cappedLine.charAt(clustaledLine.getPositionInOriginatingSequence())))
			{
				clustaledLine.incrementPositionInOriginalSequence();
			}										
		}
			
		return rSeqResidue;

	}

	private static void addChains(PdbChain[] chains,
							char[][] alignmentChars,
							int x,
							List clustalLines,
							ReferenceSequenceResidue rSeqResidue,
							Alignment a,
							String cappedLine) throws Exception
	{			
		
		for ( int y=0; y < chains.length; y++ )
		{
			ClustalLine clustalLine = ReferenceSequenceUtils.getClustalLine(clustalLines, chains[y].getChainChar());
			PdbResidue pdbResidue = null;
			
			if ( MapResiduesToIndex.isValidResidueChar(alignmentChars[clustalLine.getClustalIndex()][x]) )
			{
			
				while ( pdbResidue == null ) 
				{
					if ( clustalLine.getPositionInOriginatingSequence()> chains[y].getHighestPdbResiduePosition() ) 
						throw new Exception("Error!  Residue " + clustalLine.getPositionInOriginatingSequence()
										+ " out of bounds for " + 
											a.getAligmentID() );
					
					pdbResidue = chains[y].getPdbResidueByPdbPosition(clustalLine.getPositionInOriginatingSequence());
					clustalLine.incrementPositionInOriginalSequence();
				}
				
				if ( pdbResidue.getPdbChar() != alignmentChars[clustalLine.getClustalIndex()][x] )
					throw new Exception("Error!  Unequal chars for " + 
										a.getAligmentID() + " " + pdbResidue.getPdbChar()	
										+ " vs " + alignmentChars[clustalLine.getClustalIndex()][x] );
			}
			
			
			if ( rSeqResidue != null ) 
			{
				if ( pdbResidue != null ) 
				{
					rSeqResidue.addLinkedPdbResidue( pdbResidue );
				}
				else
				{
					pdbResidue = new PdbResidue(PdbResidue.BLANK_CHAR,
												 PdbResidue.BLANK_POSITION,
												  chains[y]);
					rSeqResidue.addLinkedPdbResidue(pdbResidue);	
				}
			}					
		}	

	}

	private static ClustalLine getAlignmentLine( List clustalList ) throws Exception
	{
		ClustalLine alignmentLine = null;
		
		for ( Iterator i = clustalList.iterator(); 
				i.hasNext();
			 )
		{
			ClustalLine cLine = (ClustalLine) i.next();
			
			if ( cLine.isAlignmentSequence() )
			{
				if ( alignmentLine!= null ) 
					throw new Exception("Error!  Dual indexes ");
					
				alignmentLine= cLine;
			}
		}
		
		return alignmentLine;	
	}

	private static ClustalLine getClustalLine( List clustalList, char chainChar ) throws Exception
	{
		ClustalLine alignmentLine = null;
		
		if ( chainChar == ' ' ) 
			chainChar = CHAIN_SPACE_CHAR;
			
		for ( Iterator i = clustalList.iterator(); 
				i.hasNext();
			 )
		{
			ClustalLine cLine = (ClustalLine) i.next();
			
			if ( ! cLine.isAlignmentSequence() && cLine.getLineName().charAt(0) == chainChar)
			{
				if ( alignmentLine!= null ) 
					throw new Exception("Error!  Dual indexes ");
					
				alignmentLine= cLine;
			}
		}
		
		
		return alignmentLine;
	}

	public static ReferenceSequence getReferenceSequence( Alignment a,
													 AlignmentLine aLine,
													 ClustalAlignment clustalAlignment,
													 BestAlignmentResultsFileLine rFileLine,
													 PdbFileWrapper pdbWrapper )
													 throws Exception
	{
		ReferenceSequence rSeq = new ReferenceSequence();
		char[][] alignmentChars = clustalAlignment.getAlignedSequences();
		
		List clustalLineList = ReferenceSequenceUtils.getClustalLines(  a, clustalAlignment, rFileLine);
		ClustalLine referenceLine = ReferenceSequenceUtils.getAlignmentLine( clustalLineList);
		PdbChain[] chains = ReferenceSequenceUtils.getPdbChains(a, clustalAlignment, rFileLine, pdbWrapper);
		
		String cappedLine = aLine.getSequence().toUpperCase();
		
		
		while( ! MapResiduesToIndex.isValidResidueChar( 
								cappedLine.charAt( referenceLine.getPositionInOriginatingSequence()) ))
			referenceLine.incrementPositionInOriginalSequence(); 
			
		for ( int x=0; x< clustalAlignment.getSequenceLength(); x++ )
		{
			ReferenceSequenceResidue rSeqResidue = ReferenceSequenceUtils.addFirstLineResidue(a, alignmentChars,
										 referenceLine,x, cappedLine);
			
			ReferenceSequenceUtils.addChains(chains, alignmentChars, x, clustalLineList, rSeqResidue, a, cappedLine);			
						
			if ( rSeqResidue != null ) 
				rSeq.addReferenceSequenceResidue(rSeqResidue);
		}
		
		
		return rSeq;	
	}					
	
	/**  The key is the 4 char pdbId
	 *    The value is a parsed and loaded pdb file, or null if that pdb file could not be found or parsed.
	 *    
	 *    if ftpNewPdbs is true, and the PDB file can't be found, the program will attempt to download it.
	 */
	public static HashMap getPdbWrappers( Alignment a, boolean ftpNewPdbs ) throws Exception
	{
		HashMap map = new HashMap();
		
		HashSet<String> pdbsToLoad = new HashSet<String>();
		
		
		for ( int x=0; x< a.getAnnotationParsers().length; x++ )
		{
			PFamPdbAnnotationParser parser = a.getAnnotationParsers()[x];
			pdbsToLoad.add( parser.getFourCharId() );
		}
		
		// try and download them if they are not there
		if ( ftpNewPdbs)
			for ( String fourCharID : pdbsToLoad )
			{
				try
				{
					PdbDownload.getOrDownloadFile(fourCharID);
				}
				catch(Exception e)
				{
					System.out.println("WARNING:  Could not download " + fourCharID);
				}
			}
		
		for ( Iterator i = pdbsToLoad.iterator();
				i.hasNext(); ) 
		{
			String fourCharId = i.next().toString();
			
			try
			{
				PdbFileWrapper fileWrapper = new PdbFileWrapper( fourCharId );
				map.put( fourCharId, fileWrapper );
			}
			catch(Exception e) 
			{
				System.out.println("------------Failed to parse " + fourCharId  );
				//e.printStackTrace();
				System.out.println("---------------------------------------------\n\n\n\n");
			}
		}
		
		return map;	
	}
	
	
	public static AlignmentPdbMatch getBestAlignmentPdbMatch(Alignment a, PdbFileWrapper fileWrapper ) throws Exception
	{
		AlignmentPdbMatch returnMatch = null;
		List list = new ArrayList();
		
		int bestLengthSoFar = -1000;
		
		for ( int x=0; x< a.getAnnotationParsers().length; x++ )
		{
			if ( a.getAnnotationParsers()[x].getFourCharId().equals(fileWrapper.getFourCharId()))
			{
			
				String pdbStringFragment = SequenceUtils.getPdbStringFragment( a.getAnnotationParsers()[x],
								 fileWrapper);
								 
				for ( Iterator i2 = a.getAlignmentLines().iterator();
							i2.hasNext(); ) 
				{
					AlignmentLine aLine = (AlignmentLine) i2.next();
					String aLineUngappedSequence = aLine.getUngappedSequence();
					int alignmentLineLength =  aLineUngappedSequence.length();
					int smallestLength = pdbStringFragment.length();
		
					if ( smallestLength > bestLengthSoFar ) 
					{
						String fasta = getFastaString(aLine, aLineUngappedSequence,
										   a.getAnnotationParsers()[x].getFourCharId() + "_" + a.getAnnotationParsers()[x].getChainChar(), 
										   	pdbStringFragment );
						ClustalAlignment clustalAlignment = ClustalWrapper.getClustalAlignment( fasta);					
						float pairwiseIdentity = clustalAlignment.getPairwiseIdentity(0,1);
						
						if ( pairwiseIdentity >= IDENTITY_CUTOFF )
						{
							System.out.println("\t\tmatched at " + pairwiseIdentity + " with length " + smallestLength );
							
							list.add( new AlignmentPdbMatch(a.getAnnotationParsers()[x], pairwiseIdentity, 
															 smallestLength, aLine));
							bestLengthSoFar = smallestLength;
						}									 		
					}
				}
			}
		}
		
		if ( list.size() == 0 ) 
			return null;
	
		if ( list.size() == 1 ) 
			return (AlignmentPdbMatch) list.get(0);
		
		Collections.sort( list );
		
		return (AlignmentPdbMatch) list.get(0);
			
	}
	
	private static List getMarkerList( Alignment a,
														   HashMap pdbWrapperMap ) throws Exception
	{
		List list = new ArrayList();
		int longestPerfect = -1;
		
		for ( int x=0; x< a.getAnnotationParsers().length; x++)
		{
			PFamPdbAnnotationParser parser= a.getAnnotationParsers()[x];
			PdbFileWrapper pdbWrapper = (PdbFileWrapper) pdbWrapperMap.get( parser.getFourCharId() );
			
			if ( pdbWrapper != null )
			{
				String pdbStringFragment = SequenceUtils.getPdbStringFragment( parser, pdbWrapper);
				
				if ( pdbStringFragment != null )
				{
					System.out.println("\tGot " + parser.getFourCharId() + "_" + parser.getChainChar() );
				
					for ( Iterator i2 = a.getAlignmentLines().iterator();
							i2.hasNext(); ) 
					{
						AlignmentLine aLine = (AlignmentLine) i2.next();
						int smallestLength = pdbStringFragment.length();
						String aLineUngappedSequence = aLine.getUngappedSequence();
						int alignmentLineLength =  aLineUngappedSequence.length();
					
						if ( alignmentLineLength < smallestLength ) 
							smallestLength = alignmentLineLength;
							
						if ( longestPerfect < smallestLength  )
						{
							String fasta = getFastaString(aLine, aLineUngappedSequence,
										 parser.getFourCharId() + "_" + parser.getChainChar(), pdbStringFragment );
							ClustalAlignment clustalAlignment = ClustalWrapper.getClustalAlignment( fasta);					
							float pairwiseIdentity = clustalAlignment.getPairwiseIdentity(0,1);
					
							if ( pairwiseIdentity > IDENTITY_CUTOFF ) 
							{
								System.out.println("\t\tmatched to " + pairwiseIdentity + " with length " + smallestLength );
						
								list.add( new AlignmentPdbMatch(parser,
															  pairwiseIdentity, 
															  smallestLength, 
															  aLine));	
															  
								if ( pairwiseIdentity == 100.0f ) 
								{
									longestPerfect = smallestLength;
									System.out.println( "\t\tlongest perfect set to " + longestPerfect );
								}
							}
						}
					}
				}
			}
		}	
		
		return list;	
	}
	
	
	private static boolean chainsAreOk(Alignment a, PdbFileWrapper wrapper) throws Exception
	{
		List list = getAllChains(a, wrapper.getFourCharId());
		
		if ( list.size() == 0 )
			throw new Exception("Logic error.  Could not find any chains for " + wrapper.getFourCharId());
			
		if ( list.size() == 1 ) 
			return true;
			
		for ( int x=0; x < list.size(); x++ ) 
		{
			PFamPdbAnnotationParser parserX = (PFamPdbAnnotationParser) list.get(x);
			String partialSequenceX = SequenceUtils.getPdbStringFragment( parserX, wrapper);	
			if ( partialSequenceX == null || partialSequenceX.trim().length() == 0 ) 
				throw new Exception("Logic error");
			
			for ( int y=x+1; y < list.size(); y++ ) 
			{
					PFamPdbAnnotationParser parserY = (PFamPdbAnnotationParser) list.get(y);
					String partialSequenceY = SequenceUtils.getPdbStringFragment( parserY, wrapper);	
					if ( partialSequenceY == null || partialSequenceY.trim().length() == 0 ) 
						throw new Exception("Logic error");
					
					String fastaString = getFasta(partialSequenceX, partialSequenceY);
					ClustalAlignment cAlignment	= ClustalWrapper.getClustalAlignment(fastaString);
					if ( cAlignment.getPairwiseIdentity(0,1) <  IDENTITY_CUTOFF )
						return false;			
			}
		}
			
		return true;
	}
	
	private static String getFasta( String string1, String string2 ) throws Exception
	{
		StringBuffer buff = new StringBuffer();
		
		buff.append(">1\n");
		buff.append( string1 + "\n" );
		buff.append(">2\n");
		buff.append( string2 + "\n" );
		
		return buff.toString();	
	}
	
	private static List getAllChains( Alignment a, String pdbId ) throws Exception
	{
		List list =new ArrayList();
		
		for ( int x=0; x< a.getAnnotationParsers().length; x++) 
		{
			PFamPdbAnnotationParser parser = a.getAnnotationParsers()[x];
			
			if ( parser.getFourCharId().equals(pdbId) )
				list.add(parser);
		}
		
		return list;	
	}
	
	public static AlignmentPdbMatch getBestMarker( Alignment a ) throws Exception
	{
		HashMap map = getPdbWrappers(a, true);
		List markerList  = getMarkerList(a, map);
		
		Collections.sort( markerList );
		
		for ( int x=0; x< markerList.size(); x++ ) 
		{
			AlignmentPdbMatch marker = (AlignmentPdbMatch) markerList.get(x);
			PdbFileWrapper wrapper = (PdbFileWrapper) map.get( marker.getParser().getFourCharId());
			
			if ( wrapper == null ) 
				throw new Exception("Logic error");
				
			if ( chainsAreOk(a, wrapper) ) 
				return marker;
			else
				System.out.println("Bad chains for " + marker.toString() );
		}
			
		return null;
	}
	
	/** aLine and aLineUngappedSequence are separte parameters for performance 
	 */
	public static String getFastaString( AlignmentLine aLine, 
																String aLineUngappedSequence,
																String pdbId,
																String pdbStringFragment ) throws Exception
	{
		StringBuffer buff = new StringBuffer();
		buff.append( ">");
		buff.append( aLine.getIdentifier());
		buff.append("\n");
		buff.append( aLineUngappedSequence);
		buff.append("\n");
		buff.append(">");
		buff.append( pdbId );
		buff.append("\n");
		buff.append( pdbStringFragment + "\n" );
		
		return buff.toString();
				
	}
								 
}