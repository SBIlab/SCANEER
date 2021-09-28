package covariance.utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;

import covariance.algorithms.ConservationSum;
import covariance.algorithms.ELSCCovariance;
import covariance.algorithms.JavaSCA;
import covariance.algorithms.McBASCCovariance;
import covariance.algorithms.OmesCovariance;
import covariance.algorithms.RandomScore;
import covariance.datacontainers.Alignment;
import covariance.datacontainers.AlignmentLine;
import covariance.datacontainers.PdbFileWrapper;
import covariance.datacontainers.PdbResidue;
import covariance.datacontainers.ReferenceSequence;
import covariance.datacontainers.ReferenceSequenceResidue;

/**  This class is just sample code.  Hopefully enough to get you started
 */
public class RunOne
{
	/** Don't use a column in the alignment if that column has gaps greater than MAX_RATIO_GAPS
	 */
	public static final float MAX_RATIO_GAPS=0.5f;

	public static void writeFileHeader( BufferedWriter writer ) throws Exception
	{
		writer.write( "alignmentI\t" );
		writer.write( "alignmentJ\t" );
		writer.write( "pdbI\t" );
		writer.write( "pdbJ\t" );
		writer.write( "cBcBDistance\t");
		writer.write( "mcBascScore\t" );
		writer.write( "omesScore\t" );
		writer.write( "javaScaScore\t" );
		writer.write( "elscScore\t" );
		writer.write( "randomScore\t" );
		writer.write( "conservationSumScore\n" );
	}
	
	public static void main(String[] args) throws Exception
	{
		// you can change these to point to files anywhere on your hard drive
		File pdbFile = new File( ConfigReader.getHomeDirectory() + File.separator 
										+  "2sns.pdb");
		File alignmentFile = new File( ConfigReader.getHomeDirectory() + 
										File.separator + "pnase.txt");
										
		
		// the line in the alignment that matches the chain in the crystal structure
		// the identifier given here must match the first token in the alignment line
		// the method findBestLine() in this class can be used to find a good line in the alignment
		String alignmentLineId = "Q99VJ0/87-221";
		
		// the chain in the crystal strucutre that matches the 
		char pdbChainChar = ' ';
		
		// the pdb residue in the pdb chain in which the alignment begins
		int pdbStartPos = 8;
		
		// the pdb residue in the pdb chain in which the alignment ends
		int pdbEndPos = 141;
		
		// remove all sequences that are greater than 90% redundant
		// this uses CLUSALW so can be slow
		// set to <0 to turn the filter off
		float filter = 90;

		// these are just used to create the results file filepath
		String familyName = "PNASE";
		File resultsFile = new File( ConfigReader.getHomeDirectory() + File.separator + 
							familyName +  "_results.txt" );
							
		BufferedWriter writer = new BufferedWriter( new FileWriter( resultsFile ));
		writeFileHeader(writer);
		
		// this pdbFileWrapper is pretty fragile
		// for example, it throws if there are any duplicate atoms
		// or if the pdb file is in some weird format
	    // this may cause you problems, depending on what you are doing
		PdbFileWrapper pdbFileWrapper = new PdbFileWrapper( pdbFile);
		
		Alignment a = new Alignment( familyName, alignmentFile, false);
		AlignmentLine aLine = a.getAnAlignmentLine(alignmentLineId);
		
		if ( aLine == null ) 
			throw new Exception( "Could not find " + alignmentLineId + " in the alignment " );
									
		RefSequenceWrapper refSeqWrapper = 
			new RefSequenceWrapper(pdbFileWrapper,
							pdbChainChar, 
							pdbStartPos,
							pdbEndPos, 
							a,
							aLine,
							filter
							);
							
		runRefSequence(refSeqWrapper, writer);
	}
	
	public static void runRefSequence( RefSequenceWrapper wrapper, BufferedWriter writer ) throws Exception
	{
		ReferenceSequence rSeq = wrapper.getReferenceSequence();
		Alignment a = wrapper.getFilteredAlignment();
		char pdbChar = wrapper.getPdbChain().getChainChar();
		
		McBASCCovariance mcBasc = new McBASCCovariance(a);
		OmesCovariance omes = new OmesCovariance(a);
		JavaSCA jSca = new JavaSCA(a);
		ELSCCovariance elsc = new ELSCCovariance(a);
		RandomScore random = new RandomScore();
		ConservationSum cSum = new ConservationSum(a);
		
		for ( int x=0; x< rSeq.getReferenceSequenceResidues().size(); x++ ) 
		{
			ReferenceSequenceResidue xResidue = 
				(ReferenceSequenceResidue) rSeq.getReferenceSequenceResidues().get(x);
			
			int xAlignmentPos = xResidue.getAlignmentPosition();
				
			if ( xResidue.allMatches() && a.getRatioValid(xAlignmentPos) > MAX_RATIO_GAPS  )
			{
				PdbResidue xPdbResidue = xResidue.getLinkedPdbResidue(pdbChar);
					
				System.out.println( "Crunching i pdb position " +  xPdbResidue.getPdbPosition());
				
				for ( int y= x+1; y < rSeq.getReferenceSequenceResidues().size(); y++ ) 
				{
					
					ReferenceSequenceResidue yResidue = 
						(ReferenceSequenceResidue) rSeq.getReferenceSequenceResidues().get(y);
						
					int yAlignmentPosition = yResidue.getAlignmentPosition();
						
					if ( yResidue.allMatches() && a.getRatioValid(yAlignmentPosition) > MAX_RATIO_GAPS ) 
					{
						PdbResidue yPdbResidue = yResidue.getLinkedPdbResidue(pdbChar);
						
						writer.write( xAlignmentPos + "\t");					
						writer.write( yAlignmentPosition + "\t");
						writer.write( xPdbResidue.getPdbPosition() + "\t" );
						writer.write( yPdbResidue.getPdbPosition() + "\t" );
						
						// write the min distance between all chains
						writer.write( xResidue.getMinCbDistance(yResidue) + "\t" );	
						writer.write( mcBasc.getScore(a, xAlignmentPos , yAlignmentPosition) + "\t" );
						writer.write( omes.getScore(a, xAlignmentPos , yAlignmentPosition) + "\t" );
						writer.write( jSca.getScore(a, xAlignmentPos , yAlignmentPosition) + "\t" );
						writer.write( elsc.getScore(a,xAlignmentPos ,yAlignmentPosition) + "\t" );
						writer.write( random.getScore(a,xAlignmentPos ,yAlignmentPosition) + "\t" );
						writer.write( cSum.getScore(a,xAlignmentPos ,yAlignmentPosition) + "\n" );
					}
					
					writer.flush();
				}
				
			}
		}
		
		writer.flush();  writer.close();
		
	}
	
	/**  Just some sample code which shows how to find the best matching alignmentLine 
	 *   between an alignment and a strucutre.
	 *   
	 *   This code is very inefficient and can hence be very slow.
	 * 
	 *   It's also more than a little ugly, but it works.
	 *   
	 *   To find the best PDB structure for a given PFAM alignment, use 
	 *   scripts.ScanForBestAlignment
	 */
	private static void findBestLine() throws Exception
	{
		File pdbFile = new File( ConfigReader.getHomeDirectory() + File.separator 
										+  "2sns.pdb");
		File alignmentFile = new File( ConfigReader.getHomeDirectory() + 
										File.separator + "pnase.txt");
										
		PdbFileWrapper wrapper = new PdbFileWrapper( pdbFile );
										
		// the chain in the crystal strucutre that matches the 
		char pdbChainChar = ' ';
		
		// the pdb residue in the pdb chain in which the alignment begins
		int pdbStartPos = 8;
		
		// the pdb residue in the pdb chain in which the alignment ends
		int pdbEndPos = 141;
										
		Alignment a = new Alignment( "anAlignment", alignmentFile, false);
	
		// a really ugly hack to make this code think your alignment came from Pfam!
		// if you want to check more than just a single chain, add more lines to the list!
		List list = new ArrayList();
		list.add( "#=GF DR   PDB; " + wrapper.getFourCharId() + " " + pdbChainChar
										+ "; " + pdbStartPos +"; " +pdbEndPos+";");
		a.setPdbIds(list);
		
		System.out.println(" Best line id is " + 
		ReferenceSequenceUtils.getBestAlignmentPdbMatch(a, wrapper ).getAlignmentLine().getIdentifier());
	}
}
