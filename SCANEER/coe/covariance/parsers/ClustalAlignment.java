package covariance.parsers;

import java.io.BufferedReader;
import java.io.BufferedWriter;

import java.io.FileWriter;
import java.io.File;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.StringTokenizer;

import covariance.datacontainers.SequenceFeature;
import covariance.utils.MapResiduesToIndex;
import covariance.utils.SequenceUtils;

public class ClustalAlignment
{
	private String originalAlignment;
	private int numSequences;
	private int sequenceLength;
	
	private char[][] alignedSequences = null;
	private int [][] charCount = null;
	private String[] sequenceNames = null;
	private List sequenceFeatures = new ArrayList();
	
	public int getNumSequences()
	{
		return this.numSequences;	
	}
	
	public int getSequenceLength()
	{
		return this.sequenceLength;
	}
	
	public String getSequence(int sequenceNum ) 
	{
		StringBuffer buff = new StringBuffer();
		
		for ( int x=0; x < alignedSequences[sequenceNum].length; x++ ) 
			buff.append( alignedSequences[sequenceNum][x] );
		
		return buff.toString();
	}
	
	public String getUngappedSequence( int sequenceNum ) 
	{
		StringBuffer buff = new StringBuffer();
		
		for ( int x=0; x < alignedSequences[sequenceNum].length; x++ ) 
			if ( alignedSequences[sequenceNum][x] != '-' ) 
				buff.append( alignedSequences[sequenceNum][x] );
		
		return buff.toString();
	}
	
	public char[][] getAlignedSequences()
	{
		return this.alignedSequences;
	}
	
	public String getSequenceName(int sequenceNum)
	{
		return this.sequenceNames[sequenceNum];
	}
	
	public String[] getSequenceNames()
	{
		return this.sequenceNames;
	}
	
	public ClustalAlignment( File alignmentFile) throws Exception
	{
		this.originalAlignment = SequenceUtils.fileToString(alignmentFile);
		
		doInitialPass();
		doSecondPass();
		fillCountArray();
	}
	
	
	private void fillCountArray()
	{
		charCount = new int[numSequences][sequenceLength];
		
		for ( int x=0; x< numSequences; x++ ) 
		{
			if ( alignedSequences[x][0] == '-' ) 
				charCount[x][0] = 0;
			else
				charCount[x][0] = 1;
		}
		
		for ( int x=0; x<numSequences; x++ ) 
			for ( int y=1; y<sequenceLength; y++) 
			{
				if ( alignedSequences[x][y] == '-' ) 
					charCount[x][y] = charCount[x][y-1];
				else
					charCount[x][y] = charCount[x][y-1] + 1;
			}
	}
	
	private void doSecondPass() throws Exception
	{
		alignedSequences = new char[ numSequences][sequenceLength];
		sequenceNames = new String[ numSequences ];
		
		BufferedReader reader = new BufferedReader( new StringReader( originalAlignment));
		
		reader.readLine();  reader.readLine();  reader.readLine();
		
		int charRead =0;
		int linesRead = 3;
		
		boolean keepGoing = true;
		String line = null;
		int baseCharRead = 0;
		
		while ( keepGoing ) 
		{	
			baseCharRead = charRead;
			
			for (int x=0; x < numSequences; x++ ) 
			{
				line = reader.readLine();
				linesRead++;
				
				StringTokenizer sToken = new StringTokenizer( line );
				String seqName = sToken.nextToken().trim();
				
				if ( sequenceNames[x] == null ) 
				{
					sequenceNames[x] = seqName;
					//System.out.println("Assigning " + seqName );
				}
				else if ( ! sequenceNames[x].equals( seqName ) ) 
				{
					throw new Exception("Error in parsing line " + linesRead + "." + 
						"Expecting sequence " + sequenceNames[x] + " but found " + seqName );
				}
				
				String charString = sToken.nextToken().trim();
				
				for ( int y = 0; y < charString.length(); y++ ) 
				{
					alignedSequences[x][y+baseCharRead] = charString.charAt(y);
					
					if ( x == 0 ) 
						charRead++;
				}
			}
			
			line = reader.readLine();
			linesRead++;
			line = reader.readLine();
			linesRead++;
			
			if ( line == null )
			{
				keepGoing = false;
				
				if ( charRead != sequenceLength )
					throw new Exception("Error in parsing line " + 
								".  Should have read " + sequenceLength + "  but instead read " + charRead );
			}
		}
	}
	
	private int countNumCharsInSequence(String inString)
	{
		StringTokenizer sToken  = new StringTokenizer(inString);
		
		sToken.nextToken();
		
		return sToken.nextToken().trim().length();
	}

	
	private void doInitialPass() throws Exception
	{
		BufferedReader reader = new BufferedReader( new StringReader( originalAlignment ));
		
		reader.readLine();  reader.readLine();  reader.readLine();
		
		numSequences = 0;
		sequenceLength = 0;
		
		List sequenceList = new ArrayList();
		
		String line = reader.readLine();
		String lastLine = null;
		
		while ( line != null && line.trim().length() != 0 && isSequenceString( line ) ) 
		{
			numSequences++;
			lastLine = line;
			line = reader.readLine();
		}
		
		if ( lastLine == null ) 
			throw new Exception("Error!  No sequence information found " );
		
		sequenceLength+= countNumCharsInSequence( lastLine );
			
		while( true )
		{
			line = reader.readLine();	
			
			if ( line == null  ) 
				return;
			
			line = reader.readLine();
			
			sequenceLength += countNumCharsInSequence( line );
			
			// read lines up to the next blank line
			for ( int x=0; x < numSequences; x++ ) 
			{
				reader.readLine();
			}
		}
	}
	
	private boolean isSequenceString( String inString ) 
	{
		if ( inString == null ) 
			return false;
		
		inString = inString.trim();
		
		if ( inString.length() ==0 )
			return false;
		
		for ( int x=0; x< inString.length(); x++) 
		{
			char c = inString.charAt(x);
			
			if(  c!=' ' && c!='\t' && c!='\n' && c!= '.' && c!='*' && c != ':' ) 
			   return true;
		}
		
		return false;
	}

	public float getMinimumPairwiseIdentity()
	{
		float min = 2000f;
		
		for ( int x=0; x< this.getNumSequences(); x++ ) 
			for ( int y=x+1; y < this.getNumSequences(); y++ ) 
			{
				float score = getPairwiseIdentity(x,y);
				
				if ( score < min)
					min= score;
			}
		
		return min;
	}

	public float getPairwiseIdentity( int sequence1, int sequence2 ) 
	{
		int length1 = 0;
		int length2 = 0;
		int numIdentical=0;
		
		for ( int x=0; x< sequenceLength; x++ ) 
		{
			char c1 = alignedSequences[sequence1][x];
			char c2 = alignedSequences[sequence2][x];
			
			if ( c1 != '-' ) 
				length1++;
			
			if ( c2 != '-' ) 
				length2++;
			
			if ( c1 == c2 )
			{
				if ( c1 != '-' || c2 != '-' ) 
					numIdentical++;
			}
		}
		
		return 100 * ((float)numIdentical)/Math.min( length1, length2 );
	}
	
	public void writeHtmlFile( File outFile ) throws Exception
	{
		BufferedWriter writer = new BufferedWriter( new FileWriter( outFile ));
		
		writer.write("<html><body><pre>");
		
		writer.write( getHTMLAlignment() );	
		
		writer.write("</pre></body></html>");		
		
		
		writer.flush();  writer.close();
		
		
		
	}
	
	/**  Will add the first and last occurence of sequenceToAdd.
	 * 
	 *   todo:  Fix so that it will add all of them 
	 */
	private void addAFeature( String sequenceToAdd,
								int sequenceNum,
								String color ) throws Exception
	{
		sequenceToAdd = SequenceUtils.stripSpaces(sequenceToAdd);
		
		String sequence = getSequence(sequenceNum);
		
		int startPosition = sequence.indexOf(sequenceToAdd);
		
		if ( startPosition == -1 ) 
			throw new Exception("Error!  Could not find " + sequenceToAdd );
			
		sequenceFeatures.add(new SequenceFeature(startPosition, startPosition + sequenceToAdd.length()-1,
						color, sequenceNum));	
						
		int endPosition = sequence.lastIndexOf(sequenceToAdd);
		
		if ( startPosition == endPosition ) 
			return;
			
		sequenceFeatures.add(new SequenceFeature(endPosition, endPosition+ sequenceToAdd.length()-1,
						color, sequenceNum));	
		
	}
	
	private String getFontTag( String oldFontString, String newFontString ) 
	{
		if ( oldFontString == null ) 
		{
			if ( newFontString == null ) 
				return "";
			
			return "<font " + newFontString + ">";
		}
		
		// if we are here, oldFoldString != null
		
		if ( newFontString == null ) 
			return "</font>";
		
		if( newFontString.equals( oldFontString ) ) 
			return "";
		else
			return "</font><font " + newFontString + ">";
	}
	
	
	
	
		/**  Returns the number of sequence characters added
	 */
	private int appendSequenceHtmlToBuffer( StringBuffer buff, String inLine, int sequenceNum, int startPosition ) 
		throws Exception
	{	
		StringTokenizer sToken = new StringTokenizer( inLine );
		String id = sToken.nextToken();
		buff.append( id );
		
		inLine = inLine.substring( id.length() );
		
		boolean keepGoing = true;
		int position = 0;
		
		while ( keepGoing ) 
		{
			if ( inLine.charAt( position ) == ' ' ) 
			{
				buff.append( ' ' );
				position++;
			}
			else
			{
				keepGoing=false;
			}
		}
		
		inLine = inLine.substring( position ).trim();
		position = 0;
		
		String fontString = null;
		String oldFontString = null;
		while( position < inLine.length() ) 
		{
			if ( inLine.charAt( position ) != alignedSequences[sequenceNum][startPosition + position] ) 
				throw new Exception("Error!  Expecting " + alignedSequences[sequenceNum][startPosition + position] + 
									" but got " + inLine.charAt( position ) );
			
			fontString = getFontString(sequenceNum, position + startPosition);
			
			buff.append( getFontTag( oldFontString, fontString ));
			
			buff.append( inLine.charAt( position ));
			
			position++;
			
			oldFontString = fontString;
		}
		
		if ( fontString != null ) 
			buff.append("</font>");
		
		buff.append("\t" + charCount[sequenceNum][startPosition + position -1] );
		
		buff.append("\n");
		
		return position;
	}

	private boolean shouldBeGreyed( int sequenceNum, int position ) 
	{
		char c = alignedSequences[sequenceNum][position];
	
		if ( c == '-' ) 
			return false;
	
		int numIdentical = 0;
		int numValid = 0;
		
		for ( int x=0; x< numSequences; x++ ) 
		{
			char targetChar = alignedSequences[x][position];
			if ( MapResiduesToIndex.isValidResidueChar(targetChar) ) 
				numValid++;
		}
		
		if ( numValid < 2 ) 
			return false;
			
		if ( numValid == 2 ) 
		{
			for ( int x=0; x< numSequences; x++ ) 
				if ( x != sequenceNum && 
						MapResiduesToIndex.isValidResidueChar(alignedSequences[x][position] )) 
				{
					if ( c == alignedSequences[x][position] )
						return true;
					else
						return false;
				}	
		}
		
		for ( int x=0; x< numSequences; x++ ) 
			if ( c == alignedSequences[x][position] ) 
				numIdentical++;
		
		if ( ((float)numIdentical) >= ((float)numValid / 2 ) ) 
			return true;
			
		return false;
	}

	private String getFeatureColor( int sequenceNum, int position ) 
	{
		for ( Iterator i = sequenceFeatures.iterator();
					i.hasNext(); ) 
		{
			SequenceFeature sFeature = (SequenceFeature) i.next();
			
			if ( sFeature.getSequenceNum() == sequenceNum &&
					position >= sFeature.getStartInclusive() &&
					  position <= sFeature.getEndInclusive() )
					  return sFeature.getColor();
		}
		
		return null;	
	}

	private String getFontString( int sequenceNum, int position )
	{
		boolean shouldBeGreyed = shouldBeGreyed( sequenceNum, position );
		String featureColor = getFeatureColor(sequenceNum, position);
		
		if ( ! shouldBeGreyed ) 
		{
			if ( featureColor == null ) 
				return null;
			else
				return "color=" + featureColor;
		}
		
		if ( featureColor == null ) 
			return "style=\"background:gray\"";
			
		
		return "style=\"background:gray\" color=" + featureColor;
	}
	
	
	public String getHTMLAlignment() throws Exception
	{
		StringBuffer buff = new StringBuffer();
		
		BufferedReader reader = new BufferedReader( new StringReader( originalAlignment ));
		
		for ( int x=0; x<3;x++) 
			buff.append( reader.readLine() + "\n" );
		
		int linesRead = 3;
		
		int baseCharRead = 0;
		int lastCharsRead = 0;
		boolean keepGoing = true;
		
		while ( keepGoing ) 
		{
			baseCharRead+=lastCharsRead;
			
			for (int x=0; x< numSequences; x++ )
			{
				lastCharsRead=appendSequenceHtmlToBuffer( buff, reader.readLine(), 	x, baseCharRead);
			}
			
			buff.append( reader.readLine() + "\n" );
			
			String nextLine = reader.readLine();
			
			if ( nextLine == null )
				keepGoing = false;
			else
				buff.append( nextLine + "\n" );
		}
		
		String shadedHTML = buff.toString();
		
		return shadedHTML;
		
	}
	
	/*
	public static void main(String[] args) throws Exception
	{
		
		ClustalAlignment cAlignment = new ClustalAlignment(
			new File( 
			"C:\\Documents and Settings\\Anthony\\Desktop\\VectorStuff\\sequences In pAcG2T\\M7_M8_inpAcG2t.aln"));
		
		cAlignment.addAFeature("GGA TCC AGT GGA AGA AAG CAC ATT GTA G", 0, "red");
		cAlignment.addAFeature( Translate.reverseTranscribe( "CTAATGATGATGATGATGATGCAGCTTCGGGGAGGTGTTGGG"), 
						0, "red");
		
		cAlignment.writeHtmlFile(new File("c:\\cygwin\\clustalw\\Seq2.html"));
		
		String aString = cAlignment.getUngappedSequence(0);
		
		aString = aString.substring( aString.indexOf("GGATCCAGTGGAAGAAAGCACATTGTAG"));
		
		System.out.println( Translate.getProteinSequence(aString));
	}
	*/
	
	public static void main(String[] args) throws Exception
	{
		
	}
	
	/*
	public static void main(String[] args) throws Exception
	{
		//String forwardSequence = Translate.readSequenceFromFastaFile(new File("C:\\cygwin\\clustalw\\58-MF.seq"));
		//System.out.println( Translate.reverseTranscribe(forwardSequence));
		
		ClustalAlignment cAlignment = new ClustalAlignment(
			new File("C:\\cygwin\\clustalw\\58-MF.aln"));
			
		
		cAlignment.addAFeature("GGA TCC AGT GGA AGA AAG CAC ATT GTA G", 0, "red");
		cAlignment.addAFeature(Translate.reverseTranscribe("CTAATGATGATGATGATGATGATCTTCAACTTCTCTGATTGG"), 0, "red");
		cAlignment.writeHtmlFile(
			new File("C:\\cygwin\\clustalw\\58-MF.html"));
			
	}
	*/
	
	/*
	public static void main(String[] args) throws Exception
	{
		ClustalAlignment cAlignment = new ClustalAlignment(
			new File("C:\\cygwin\\clustalw\\dros.aln"));
		cAlignment.writeHtmlFile(new File("C:\\cygwin\\clustalw\\dros.html"));
	}
	*/
	
	/*
	 *
	 *
	 * public static void main(String[] args) throws Exception
	{
		ClustalAlignment cAlignment = new ClustalAlignment(
			new File("C:\\Documents and Settings\\Anthony\\Desktop\\VectorStuff\\sequencesInTopo\\D_Full.aln"));
			
		System.out.println( cAlignment.getSequenceName(2) );
		String aSequence = cAlignment.getUngappedSequence(2);
		aSequence = aSequence.substring(aSequence.indexOf("GGATCCGAAC") );
		aSequence = aSequence.substring(0, aSequence.indexOf('N'));
		System.out.println( aSequence );
		String translatedRealSequence = Translate.getProteinSequence(aSequence);
		System.out.println( translatedRealSequence);
		
		String fullSequence = cAlignment.getUngappedSequence(0);
		fullSequence = fullSequence.substring(fullSequence.indexOf("GGATCCGAAC"));
		String translatedFullSequence = Translate.getProteinSequence(fullSequence);
		System.out.println( translatedFullSequence );
		translatedFullSequence = translatedFullSequence.substring(0, translatedFullSequence.indexOf("EVKRAWFYCKAC"));
		
		for ( int x=0; x < translatedFullSequence.length(); x++ )
		{
			if ( translatedFullSequence.charAt(x) != translatedRealSequence.charAt(x) ) 
				System.out.println( "Diff " + x + " " + translatedFullSequence.charAt(x) + " " + 
								translatedRealSequence.charAt(x)  );
		}		
	}*/
}
