package covariance.datacontainers;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.StringTokenizer;

import covariance.algorithms.ConservationSum;
import covariance.parsers.PFamPdbAnnotationParser;
import covariance.utils.MapResiduesToIndex;


public class Alignment
{
	public static final Character GAP_CHARACTER = new Character( '-');
	public static final char GAP_CHAR= GAP_CHARACTER.charValue();
	
	private String alignmentID;
	private List alignmentLines;
	private int numColumns=-1;
	private int[][] counts = null;
	private float[] frequencies = null;
	private int[] totalValid = null;
	private PFamPdbAnnotationParser pdbIds[];
	private String[] columnStrings;
	
	public String getAligmentID()
	{
		return alignmentID;
	}
	
	public boolean columnHasValidResidue( int col ) throws Exception
	{
		if ( getRatioValid(col) > 0f ) 
			return true;
			
		return false;
	}

	
	/**  See Henikoff and Henikoff, J. Mol. Biol. 1994, 243, 574-578.
	 * 
	 *   Returns the array of weights, which probably will only be used for testing
	 **/
	public synchronized double[] adjustFrequenciesByWeights() throws Exception
	{
		getFrequencies();  // initialize
		
		double[] weights = new double[ getNumSequencesInAlignment() ];
		double weightSum = 0;
		
		for ( int col=0; col < getNumColumnsInAlignment(); col++ ) 
		{
			int r =0;
			
			for ( int res=0; res < MapResiduesToIndex.NUM_VALID_RESIDUES; res++ )
					if ( counts[col][res] > 0 ) 
						r++;	
						
			System.out.println("Got r of " + r + " for col " + col);
			
			for ( int row=0; row< getNumSequencesInAlignment(); row++ ) 
			{
				String aSequence = ((AlignmentLine) getAlignmentLines().get(row)).getSequence();
				char c = aSequence.charAt(col);
				
				if ( MapResiduesToIndex.isValidResidueChar(c) ) 
				{
					weights[row] += 1f / ( r * counts[col][MapResiduesToIndex.getIndex(c)] );
				}
			}
		}
		
	
		for ( int x=0; x< weights.length; x++ ) 
			weightSum += weights[x];
		
		for ( int x=0; x< weights.length; x++ ) 
		{
			weights[x] /= weightSum;
		}
			
		this.frequencies = new float[ MapResiduesToIndex.NUM_VALID_RESIDUES];
		int totalNum =0;
			
		for ( int res = 0; res < MapResiduesToIndex.NUM_VALID_RESIDUES; res++ ) 
		{
			
			for ( int row = 0; row < getNumSequencesInAlignment(); row++ ) 
			{
				String aSequence = ((AlignmentLine) getAlignmentLines().get(row)).getSequence();
				
				for ( int col=0; col < getNumColumnsInAlignment(); col++ ) 
				{
					char c = aSequence.charAt(col);
					
					if ( MapResiduesToIndex.isValidResidueChar(c) )
					{
						totalNum++;
						this.frequencies[res] += weights[row] * getNumSequencesInAlignment();
					}			
				}
			}
				
		}
		
		for ( int x=0; x< MapResiduesToIndex.NUM_VALID_RESIDUES; x++ )
			this.frequencies[x] /= totalNum;
		
			
		return weights;
	}
	
	public float getMeanConservation( float cutoffRatio ) throws Exception
	{
		ConservationSum cSum = new ConservationSum(this);
		float sum = 0;
		int n =0;
		
		for ( int x=0; x< this.numColumns; x++ )
		{
			if ( getRatioValid(x) >= cutoffRatio )
			{
				sum+= cSum.getScore(x);
				n++;
			}
		}
		
		return sum / n;
	}
	
	
	public String getColumnAsString( int colPos ) throws Exception
	{
		return columnStrings[colPos];
	}
	
	/**   removes any column in which the most conserved residue in that column 
	 *    is less than the cutoff ( cutoff should be between 0 and 1 ) 
	 */
	public Alignment getAlignmentWithRemovedUnconservedColumns(float cutoff ) throws Exception
	{
		boolean[] tossColumns = getTossColumns(cutoff);
		List newAlignmentLines = getNewAlignmentLines(tossColumns);
		
		Alignment newAlignment = new Alignment(this.getAligmentID(), newAlignmentLines);	
		newAlignment.setPdbIds( this.getPdbIds());
		
		return newAlignment;
	}
	
	public static Alignment getAlignmentFromFasta( File file ) throws Exception
	{
		BufferedReader reader = new BufferedReader( new FileReader( file ));
		List list= new ArrayList();
		
		String nextLine = reader.readLine();
		
		while ( nextLine != null ) 
		{
			list.add( new AlignmentLine( nextLine, reader.readLine() ));
			nextLine = reader.readLine();
		}
		
		return new Alignment("fromFasta", list );
	}
	
	private List getNewAlignmentLines( boolean[] tossColumns ) 
	{
		List newAlignmentLines = new ArrayList();
		
		for ( Iterator i = this.getAlignmentLines().iterator();
				i.hasNext(); ) 
		{
			AlignmentLine aLine = (AlignmentLine) i.next();
			
			StringBuffer buff = new StringBuffer();
			String oldSequence = aLine.getSequence();
			
			for ( int x=0; x < oldSequence.length(); x++ ) 
				if ( ! tossColumns[x] ) 
					buff.append(oldSequence.charAt(x));
					
			newAlignmentLines.add( new AlignmentLine( aLine.getIdentifier(), buff.toString()));
		}
		
		return newAlignmentLines;
	}
	
	/**  Each tossColumn returns true if the most frequent residue in the column occurs
	 *   in less than the cutoff.
	 * 
	 *   For example, if the cutoff value is .2, any column in which the most conserved residue
	 *   is less than 20% of all valid residues in that column, will have a tossColumn of false.
	 * 
	 *   This is used for pertubation based covariation algorithms which require a minimal number 
	 *   of sequences in their pertubation.
	 */
	public boolean[] getTossColumns( float cutoff ) throws Exception
	{
		int[][] counts = getCounts();
		
		if ( cutoff < 0 || cutoff > 1 ) 
			throw new Exception("Illegal cutoff " + cutoff );
			
		char[] mostFrequentChars = getMostFrequentResidues(this);
		boolean[] tossColumn = new boolean[this.numColumns];
		
		for ( int x=0; x< this.numColumns; x++ ) 
		{
			if ( ! MapResiduesToIndex.isValidResidueChar(mostFrequentChars[x]) )
			{
				tossColumn[x] = true;
			}
			else
			{
				float ratio = 
					((float)counts[x][MapResiduesToIndex.getIndex(mostFrequentChars[x])]) 
													/ this.getNumSequencesInAlignment();
				
				if ( ratio > cutoff ) 
					tossColumn[x] = false;
				else
					tossColumn[x] = true;
			} 
		}
		
		return tossColumn;
	}
	
	/**  Tosses the second of any two sequences that have a percent identity greater than or equal to the 
	 *   cutoff ( cutoff should be between 0 and 100 )
	 */ 
	public Alignment getFilteredAlignment(float cutoff) throws Exception
	{
		
		if ( cutoff > 100 ) 
			throw new Exception("Error!  Illegal cutoff " + cutoff);	
		
		List newLines = new ArrayList();
		
		for ( Iterator i = this.getAlignmentLines().iterator();
				i.hasNext(); ) 
		{
			AlignmentLine aLine = (AlignmentLine) i.next();
			
			if ( okToAdd(newLines, aLine, cutoff) ) 
				newLines.add(aLine);
		}
		
		return new Alignment( this.getAligmentID(), newLines );		
	}
	
	private boolean okToAdd( List lines, AlignmentLine aLine, float cutoff ) 
	{
		for ( Iterator i = lines.iterator();
				i.hasNext(); ) 
		{
			AlignmentLine aListLine = (AlignmentLine) i.next();
			
			if ( getPairwiseIdentity( aLine, aListLine ) >= cutoff )
				return false; 
		}
		
		return true;	
	}
	
	public void writeUngappedSequencesAsFasta( File file ) throws Exception
	{
		BufferedWriter writer = new BufferedWriter( new FileWriter( file ));
		writeUngappedSequencesAsFasta(writer);
		writer.flush();  writer.close();	
	}
	
	public void writeUngappedSequencesAsFasta( BufferedWriter writer ) throws Exception
	{
		for ( Iterator i = alignmentLines.iterator();
			  i.hasNext(); ) 
		{
			AlignmentLine aLine = (AlignmentLine) i.next();
			
			writer.write( ">" + aLine.getIdentifier() + "\n" );
			writer.write( aLine.getUngappedSequence().toString() );
			writer.write( "\n");
		}
		
		writer.flush();
	}
	
	public int getNumSequencesInAlignment()
	{
		return alignmentLines.size();
	}
	
	public int getNumColumnsInAlignment()
	{
		return numColumns;
	}
	
	public List getAlignmentLines()
	{
		return this.alignmentLines;
	}
	
	public static float getAveragePairwiseIdentity( List alignmentLines ) 
	{
		float sum = 0;
		int n =0;
		
		for ( int x=0; x< alignmentLines.size(); x ++ ) 
		{
			AlignmentLine xLine = (AlignmentLine) alignmentLines.get(x);
			
			for ( int y=x+1; y < alignmentLines.size(); y++) 
			{
				AlignmentLine yLine = (AlignmentLine) alignmentLines.get(y);
				
				sum += getPairwiseIdentity(xLine, yLine);
				n++;	
			}
		}
		
		return sum / n;
	}
	
	public static float getPairwiseIdentity( AlignmentLine line1, AlignmentLine line2) 
	{
		int length1 = 0;
		int length2 = 0;
		String sequence1 = line1.getSequence().toUpperCase();
		String sequence2 = line2.getSequence().toUpperCase();
		int numIdentical=0;
		
		for ( int x=0; x< line1.getSequence().length(); x++ ) 
		{
			char c1 = sequence1.charAt(x);
			char c2 = sequence2.charAt(x);
			boolean c1Ok = MapResiduesToIndex.isValidResidueChar(c1);
			boolean c2Ok = MapResiduesToIndex.isValidResidueChar(c2);
			
			if ( c1Ok) 
				length1++;
			
			if ( c2Ok ) 
				length2++;
			
			if ( c1Ok && c2Ok && c1 == c2 )
			{
				numIdentical++;
			}
		}
		
		return 100 * ((float)numIdentical)/Math.min( length1, length2 );
	}		
	
	public boolean containsAnAlignmentLine( String id )
	{
		for ( Iterator i = getAlignmentLines().iterator();
				i.hasNext(); ) 
		{
			AlignmentLine aLine = ( AlignmentLine) i.next();
			
			if ( aLine.getIdentifier().equals(id) )
				return true;
		}
		
		return false;
		
	}
	
	public AlignmentLine getAnAlignmentLine(String id) throws Exception
	{
		for ( Iterator i = getAlignmentLines().iterator();
				i.hasNext(); ) 
		{
			AlignmentLine aLine = ( AlignmentLine) i.next();
			
			if ( aLine.getIdentifier().equals(id) )
				return aLine;
		}
		
		throw new Exception("Could not find " + id );
	}
	
	public int getTotalNumValidColumns( float cutoff ) throws Exception
	{
		int sum = 0;
		
		for ( int x=0; x< getNumColumnsInAlignment(); x++ ) 
			if ( getRatioValid(x) >= cutoff ) 
				sum++;
		
		return sum;
	}
	
	public synchronized int[] getTotalNumValidResidues() throws Exception
	{
		if ( totalValid == null ) 
		{
			int[][] counts = getCounts();
			totalValid = new int[ counts.length];
			
			for ( int x=0; x < totalValid.length; x++ ) 
				totalValid[x] = 0;
			
			for (int x=0; x< counts.length; x++) 
			{
				for ( int y=0; y < MapResiduesToIndex.NUM_VALID_RESIDUES; y++ ) 
				{
					totalValid[x] += counts[x][y];
				}
			}
		}
		
		return totalValid;
	}
	
	/**  returns an alignment line id whose ungapped sequence is a perfect subfragment of the passed in sequence,
	 *   if one exists.
	 * 
	 *   Otherwise returns null
	 */
	public AlignmentLine findPerfectSubframgent( String fragment ) throws Exception
	{
		fragment = fragment.toUpperCase();
		StringBuffer queryBuff = new StringBuffer();
		
		for ( int x=0; x< fragment.length(); x++ ) 
		{
			char c= fragment.charAt(x);
			
			if ( MapResiduesToIndex.isValidResidueChar(c))
				queryBuff.append(c);		
		}
		
		fragment = queryBuff.toString();
		
		for ( Iterator i = this.getAlignmentLines().iterator();
				i.hasNext(); ) 
		{
			AlignmentLine aLine = (AlignmentLine) i.next();
			String theSequence = aLine.getUngappedSequence().toUpperCase();
			
			if ( theSequence.indexOf(fragment) != -1 ||
					fragment.indexOf(theSequence) != -1 )
					return aLine;			
		}
			
		return null;
	}
	
	/**  The first token in the file should be an identifier
	 *   The second token should be the sequence
	 */
	
	public static List getAlignmentLinesFromFile( File file, boolean forceToUpperCase ) throws Exception
	{
		List alignmentLines = new ArrayList();
		
		BufferedReader reader = new BufferedReader( new FileReader( file ));
		
		String nextLine = reader.readLine();
		
		while ( nextLine != null && nextLine.trim().length() > 0 ) 
		{
			StringTokenizer sToken = new StringTokenizer( nextLine );
			
			String id = sToken.nextToken();
			String sequence = sToken.nextToken();
			
			if ( forceToUpperCase ) 
				sequence = sequence.toUpperCase();
			
			alignmentLines.add( new AlignmentLine(id, sequence));
			//System.out.println( sequence );
			
			nextLine = reader.readLine();
		}
		
		return alignmentLines;
	}
	
	/**  The first token in the file should be an identifier
	 *   The second token should be the sequence
	 */
	public Alignment( String alignmentID,
					   File file,
					   boolean forceToUpperCase) throws Exception
	{
		this(alignmentID, getAlignmentLinesFromFile(file, forceToUpperCase ));
	}
	
	public Alignment( String alignmentID,
					  List alignmentLines) throws Exception
	{
		this.alignmentID = alignmentID;
		this.alignmentLines = alignmentLines;
		this.numColumns = assertAllAlignmentsEqualLength();	
		
		this.columnStrings = new String[ this.numColumns];
		
		for ( int x=0; x < this.numColumns; x++ )
		{
		
			StringBuffer buff = new StringBuffer();
		
			for ( Iterator i = this.alignmentLines.iterator();
					i.hasNext(); ) 
			{
				AlignmentLine aLine = (AlignmentLine) i.next();
				buff.append( aLine.getSequence().charAt(x));				
			}
		
			columnStrings[x] = buff.toString();
		}
	}
	
	/**  returns the length of all of the alignments
	 */
	private int assertAllAlignmentsEqualLength() throws Exception
	{
		if ( this.alignmentLines.size() == 0 ) 
			return 0;
		
		int aLength = ((AlignmentLine) this.alignmentLines.iterator().next()).getSequence().length();
		
		for ( Iterator i = this.alignmentLines.iterator();
			  i.hasNext();
			  )
		{
			AlignmentLine aLine = (AlignmentLine) i.next();
			//System.out.println(aLine.getIdentifier() + " " + aLine.getSequence() );
			
			String sequence = aLine.getSequence();
			
			if ( sequence.length() != aLength )
				throw new Exception("Expecting a length of " + aLength + " but got a length of " + 
											sequence.length() );

		}
		
		return aLength;
	}
	
	
	/**  The first index goes from 0 to i, representing columns in the multiple sequence alignment
	 *   the second index goes from 0 to 19 representing types of Amino Acids. 
	 *   Use MapResiduesToIndex.getChar() to map from the numbers to single-letter amino acid codes
	 */
	public synchronized int[][] getCounts() throws Exception
	{
		if ( this.counts == null ) 
		{
			
			this.counts = new int[ numColumns][ MapResiduesToIndex.NUM_VALID_RESIDUES];
			
			for ( int x=0; x < numColumns; x++ ) 
			{	
				for ( int y =0; y < alignmentLines.size(); y++ ) 
				{
					AlignmentLine aLine = (AlignmentLine) alignmentLines.get(y);	
					String sequence = aLine.getSequence();
					
					// cycle through all 20 residues
					for (int z=0; z < MapResiduesToIndex.NUM_VALID_RESIDUES; z++ ) 
					{
						if ( sequence.charAt(x) == MapResiduesToIndex.getChar(z) ) 
							counts[x][z]++;
						
					}
				}
			}
			
		}
		
		return this.counts;
	}
	
	private String getFirstTwentyChars( String inString ) 
	{
		if ( inString.length() >= 20 ) 
			return inString.substring( 0, 20 );
		
		StringBuffer buff = new StringBuffer();
		buff.append( inString );
		
		while ( buff.length() < 20 ) 
			buff.append( " " );
		
		return buff.toString();
	}
	
	/**  Returns a number between 0 and 1 indicating the fraction of residues that are valid in 
	 *    column i.  (0 == none valid, 1 == all valid ).
	 */
	public float getRatioValid( int col) throws Exception
	{
		int[][] counts = getCounts();
		int sum = 0;
		
		for (int x=0; x < MapResiduesToIndex.NUM_VALID_RESIDUES; x++ ) 
			sum+= counts[col][x];
		
		return ((float) sum)/this.alignmentLines.size();
	}
	
	/**  Returns an array that with 20 frequencies corresponding to the relative frequency of 
	 *   each residue in the MSA
	 */
	public synchronized float[] getFrequencies() throws Exception
	{
		if ( this.frequencies == null ) 
		{
			getCounts();  // initialize
			
			this.frequencies = new float[ MapResiduesToIndex.NUM_VALID_RESIDUES];
			
			for ( int x=0; x < MapResiduesToIndex.NUM_VALID_RESIDUES; x++ ) 
				this.frequencies[x] = 0;
			
			for ( int y=0; y< MapResiduesToIndex.NUM_VALID_RESIDUES; y++ ) 
				for ( int x=0; x< counts.length; x++ ) 
					this.frequencies[y] += counts[x][y];
			
			float totalNum = 0;
			
			for ( int x=0; x < MapResiduesToIndex.NUM_VALID_RESIDUES; x++ ) 
				totalNum += this.frequencies[x];
			
			for ( int x=0; x < MapResiduesToIndex.NUM_VALID_RESIDUES; x++ ) 
				this.frequencies[x] /= totalNum;
		}
		
		return this.frequencies;
	}
	
	public void dumpSubsetOfColsAsHtml( File file, int column1Pos, char column1Char,
											int column2Pos, char column2Char ) throws Exception
	{
		if ( column1Pos > column2Pos ) 
			throw new Exception("Columns out of order");
		
		BufferedWriter writer = new BufferedWriter( new FileWriter( file ));
		writer.write("<html><body>");
		writer.write("<pre>");		
		
		for ( int x=0; x<20; x++ ) 
			writer.write("-");
			
		writer.write( "" + ( column1Pos/ 100 ));
		writer.write( "" + ( column2Pos/ 100 ));
		
		writer.write( "\n");
		
		for ( int x=0; x<20; x++ ) 
			writer.write("-");
			
		writer.write( "" + ( ( column1Pos - (column1Pos / 100) * 100 ) / 10));
		writer.write( "" + ( ( column2Pos - (column2Pos / 100) * 100 ) / 10));
			
		writer.write("\n");
		
		for ( int x=0; x<20; x++ ) 
			writer.write("-");
		
		writer.write( "" + ( column1Pos - (column1Pos/10) * 10 ));
		writer.write( "" + ( column2Pos - (column2Pos/10) * 10 ));
		
		writer.write("\n\n");
		
		for ( Iterator i=this.alignmentLines.iterator();
			  i.hasNext(); ) 
		{
			AlignmentLine aLine = (AlignmentLine) i.next();
			
			String aSequence = aLine.getSequence();
			writer.write( getFirstTwentyChars( aLine.getIdentifier() ));
	
			
			for ( int x=0; x< aSequence.length(); x++ ) 
			{
				if ( x == column1Pos || x == column2Pos  )
				{
					boolean addColor =false;
					
					if (  x == column1Pos && aSequence.charAt(x) == column1Char  )
						addColor = true;
					
					if ( x == column2Pos && aSequence.charAt(x) == column2Char ) 
						addColor = true;		  
						
					if ( addColor ) 
						writer.write("<font face==\"Courier\" color=red>");
						
					writer.write( aSequence.charAt(x));
						
					if ( addColor ) 
						writer.write( "</font>");		
				}
			}
			
			writer.write( "\n");
		}
		
		
		writer.write("</pre></body></html>");
		writer.flush();  writer.close();
	}										
	
	/**  only shows 10 residues on either side of pairs
	 */
	public void dumpSubsetAlignmentAsHTML( File file, List pairColors ) throws Exception
	{
		
	}
	
	public void dumpUngappedAlignmentAsHtml( File file, List pairColors ) throws Exception
	{
		getCounts(); // initialzie counts
		BufferedWriter writer = new BufferedWriter( new FileWriter( file ));
		writer.write("<html><body>");
		writer.write("<pre>");		
		
		dumpAlignmentHeaderWithStrippedGaps(writer);
		
		for ( Iterator i=this.alignmentLines.iterator();
			  i.hasNext(); ) 
		{
			AlignmentLine aLine = (AlignmentLine) i.next();
			
			String aSequence = aLine.getSequence();
			writer.write( getFirstTwentyChars( aLine.getIdentifier() ));
	
			
			for ( int x=0; x< aSequence.length(); x++ ) 
			{
				if ( getRatioValid(x) > 0.5 ) 
				{
				
					boolean addColor =false;
					boolean makeGrey = false;
					
					if ( MapResiduesToIndex.isValidResidueChar(aSequence.charAt(x)) &&
							100f * counts[x][MapResiduesToIndex.getIndex(aSequence.charAt(x))]
								/ getNumSequencesInAlignment() >= 40f ) 
							makeGrey = true;
					
					String colorString = PairColors.getColorStringAsRange(pairColors, aLine, x );
					
					if ( colorString != PairColors.DEFAULT_COLOR_STRING ) 
						addColor = true;
						
					if ( addColor || makeGrey ) 
						writer.write( "<font ");
						
					if ( addColor ) 
						writer.write(" color= " + colorString + " " );
						
					if ( makeGrey ) 
						writer.write( " style=\"background:gray\" ");
						
					if ( addColor || makeGrey ) 
						writer.write( ">");
					
					writer.write( aSequence.charAt(x));
						
					if ( addColor || makeGrey) 
						writer.write( "</font>");
						
				}
			}
			
			writer.write("     " + aLine.getIdentifier() );
			writer.write( "\n");
		}
		
		
		writer.write("</pre></body></html>");
		writer.flush();  writer.close();		
	}
	
	public void dumpAlignmentAsHTML( File file, List pairColors ) throws Exception
	{
		getCounts(); // initialzie counts
		BufferedWriter writer = new BufferedWriter( new FileWriter( file ));
		writer.write("<html><body>");
		writer.write("<pre>");		
		
		dumpAlignmentHeader(writer);
		
		for ( Iterator i=this.alignmentLines.iterator();
			  i.hasNext(); ) 
		{
			AlignmentLine aLine = (AlignmentLine) i.next();
			
			String aSequence = aLine.getSequence();
			writer.write( getFirstTwentyChars( aLine.getIdentifier() ));
	
			
			for ( int x=0; x< aSequence.length(); x++ ) 
			{
				boolean addColor =false;
				boolean makeGrey = false;
				
				if ( MapResiduesToIndex.isValidResidueChar(aSequence.charAt(x)) &&
						100f * counts[x][MapResiduesToIndex.getIndex(aSequence.charAt(x))]
							/ getNumSequencesInAlignment() >= 40f ) 
						makeGrey = true;
				
				String colorString = PairColors.getColorString(pairColors, aLine, x );
				
				if ( colorString != PairColors.DEFAULT_COLOR_STRING ) 
					addColor = true;
					
				if ( addColor || makeGrey ) 
					writer.write( "<font ");
					
				if ( addColor ) 
					writer.write(" color= " + colorString + " " );
					
				if ( makeGrey ) 
					writer.write( " style=\"background:gray\" ");
					
				if ( addColor || makeGrey ) 
					writer.write( ">");
				
				writer.write( aSequence.charAt(x));
					
				if ( addColor || makeGrey) 
					writer.write( "</font>");
			}
			
			writer.write( "\n");
		}
		
		
		writer.write("</pre></body></html>");
		writer.flush();  writer.close();
		
	}
	
	public void dumpAlignmentHeader( BufferedWriter writer ) throws Exception
	{
		for ( int x=0; x<20; x++ ) 
			writer.write("-");
			
		for ( int x=1; x<=getNumColumnsInAlignment(); x++ ) 
			writer.write( "" + ( x / 100 ));
		
		writer.write( "\n");
		
		for ( int x=0; x<20; x++ ) 
			writer.write("-");
			
		for ( int x=1; x<=getNumColumnsInAlignment(); x++ ) 
			writer.write( "" + ( ( x - (x / 100) * 100 ) / 10));
			
		writer.write("\n");
		
		for ( int x=0; x<20; x++ ) 
			writer.write("-");
		
		for ( int x=1; x<=getNumColumnsInAlignment(); x++ ) 
		{
			writer.write( "" + ( x - (x/10) * 10 ));
		}
		
		writer.write("\n\n");
			
	}
	
	private void dumpAlignmentHeaderWithStrippedGaps( BufferedWriter writer ) throws Exception
	{
		for ( int x=0; x<20; x++ ) 
			writer.write("-");
			
		for ( int x=1; x<=getNumColumnsInAlignment(); x++ ) 
			if ( getRatioValid(x-1) > 0.5 ) 
				writer.write( "" + ( x / 100 ));
		
		writer.write( "\n");
		
		for ( int x=0; x<20; x++ ) 
			writer.write("-");
			
		for ( int x=1; x<=getNumColumnsInAlignment(); x++ ) 
			if ( getRatioValid(x-1) > 0.5 ) 
				writer.write( "" + ( ( x - (x / 100) * 100 ) / 10));
			
		writer.write("\n");
		
		for ( int x=0; x<20; x++ ) 
			writer.write("-");
		
		for ( int x=1; x<=getNumColumnsInAlignment(); x++ ) 
			if ( getRatioValid(x-1) > 0.5 ) 
				writer.write( "" + ( x - (x/10) * 10 ));
		
		writer.write("\n\n");	
	}
	
	public void dumpAlignmnetWithHeader( File file ) throws Exception
	{
		BufferedWriter writer = new BufferedWriter( new FileWriter( file ));
		
		dumpAlignmentHeader(writer);
		dumpAlignment(writer);			
		writer.flush();  writer.close();	
	}
	
	public void dumpAlignmentAsFasta( File file ) throws Exception
	{
		BufferedWriter writer = new BufferedWriter( new FileWriter( file ));
		
		for ( Iterator i = alignmentLines.iterator();
				i.hasNext(); ) 
		{
			AlignmentLine aLine = (AlignmentLine) i.next();
			writer.write( ">" + aLine.getIdentifier() + "\n" );
			writer.write( aLine.getSequence() + "\n" );
		}
		
		writer.flush();  writer.close();
	}
	
	public void dumpAlignment( File file ) throws Exception
	{
		BufferedWriter writer = new BufferedWriter( new FileWriter( file ));
		dumpAlignment(writer);
		writer.flush();  writer.close();
	}
	
	public void dumpAlignment( BufferedWriter writer ) throws Exception
	{
		for ( Iterator i=this.alignmentLines.iterator();
			  i.hasNext(); ) 
		{
			AlignmentLine aLine = (AlignmentLine) i.next();
			writer.write( getFirstTwentyChars( aLine.getIdentifier() ));
			writer.write( aLine.getSequence() );
			writer.write( "\n");
		}
		
		writer.flush();
	}
	
	public PFamPdbAnnotationParser[] getAnnotationParsers()
	{
		return pdbIds;
	}
	
	public PFamPdbAnnotationParser getAPdbIdParser(String pdbId, char chain)
	{
		for ( int x=0; x< this.pdbIds.length; x++ )
		{
			if ( pdbIds[x].getChainChar() == chain &&
			  	 pdbIds[x].getFourCharId().equals(pdbId) ) 
			  		return pdbIds[x];
		}
		
		return null;
	}
		
	private void setPdbIds( PFamPdbAnnotationParser[] newPdbIDs )
	{
		this.pdbIds = newPdbIDs;
	}
		
	public void setPdbIds(List pdbStringList) throws Exception
	{
		this.pdbIds = new PFamPdbAnnotationParser[ pdbStringList.size() ];
		int num = 0;
		
		for ( Iterator i = pdbStringList.iterator();
				i.hasNext(); ) 
		{
			pdbIds[num] = new PFamPdbAnnotationParser( i.next().toString() );			
			num++;
		}
	}

	public PFamPdbAnnotationParser[] getPdbIds()
	{
		return pdbIds;
	}

	public static char[] getMostFrequentResidues(Alignment anAlignment) throws Exception
	{
		char[] residues = new char[ anAlignment.getNumColumnsInAlignment() ];
		int[][] counts = anAlignment.getCounts();
		
		for ( int x=0; x< anAlignment.getNumColumnsInAlignment(); x++ ) 
		{
			long max = 0;
			residues[x] = '-';
			
			for ( int y=0; y < MapResiduesToIndex.NUM_VALID_RESIDUES; y++ ) 
			{
				if ( counts[x][y] > max ) 
				{
					max = counts[x][y];
					residues[x] = MapResiduesToIndex.getChar( y );
				}
			}
		}
		
		return residues;
	}
	


}
