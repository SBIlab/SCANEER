package scripts;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.HashSet;
import java.util.StringTokenizer;

import covariance.datacontainers.Alignment;
import covariance.datacontainers.AlignmentPdbMatch;
import covariance.parsers.PfamParser;
import covariance.utils.ConfigReader;
import covariance.utils.ReferenceSequenceUtils;

public class ScanForBestAlignment 
{
	public static HashSet<String> getIncludedFamilies() throws Exception
	{
		HashSet<String> set = new HashSet<String>();
		
		File includedFamiliesFile = 
			new File(
					 ConfigReader.getHomeDirectory() + File.separator  + 
					 "includedFamilies.txt");
		
		if ( ! includedFamiliesFile.exists() )
		{
			System.out.println("WARNING:  Could not find " +
					includedFamiliesFile.getAbsolutePath() + ".  Including all families");
			return null;
		}
		
		BufferedReader reader = new BufferedReader(new FileReader(includedFamiliesFile));
		
		String nextLine = reader.readLine();
		
		while ( nextLine != null)
		{
			set.add(new StringTokenizer(nextLine).nextToken().toUpperCase());
			nextLine = reader.readLine();
		}
		
		reader.close();
		
		return set;
	}
	
	public static void main(String[] args) throws Exception
	{
		HashSet<String> includedFamilies = getIncludedFamilies();
		int alignmentNum = 0;
		PfamParser pFamParser = new PfamParser();
		
		File outFile = new File( ConfigReader.getHomeDirectory() + 
				File.separator + "AlignmentScanResults.txt" );
		System.out.println("writing to " + outFile.getAbsolutePath() );
		FileWriter writer = new FileWriter(outFile );
											
		writer.write( "id\talignmentNum\tnumColumnsInAlignment\tnumSequencesInAlignment\t" );
		writer.write("alingmentLineId\tungappedResidueLength\tpdbId\tpdbChain\tpdbStartResidue\tpdbEndResidue\tpdbLength\tpercentIdentity\tminimunLength\n");	
		writer.flush();
		
		while( ! pFamParser.isFinished()) 
		{
			alignmentNum++;
			
			Alignment a = pFamParser.getNextAlignment();
			
			if ( a != null)
			{
				System.out.println("Parsed " + alignmentNum + " alignments.  Last was " + a.getAligmentID() );
			
				if ( includedFamilies == null || includedFamilies.contains(a.getAligmentID().toUpperCase())) 
				{
					
					AlignmentPdbMatch marked = ReferenceSequenceUtils.getBestMarker(a);
					
					if ( marked != null ) 
					{
						writer.write( a.getAligmentID() + "\t" );
						writer.write( alignmentNum + "\t");
						writer.write( a.getNumColumnsInAlignment() + "\t" );
						writer.write( "" + a.getNumSequencesInAlignment() );
						writer.write( "\t" + marked.getAlignmentLine().getIdentifier() );
						writer.write( "\t" + marked.getAlignmentLine().getUngappedSequence().length()  );
						writer.write( "\t" + marked.getParser().getFourCharId() );
						writer.write( "\t" + marked.getParser().getChainChar());
						writer.write( "\t" + marked.getParser().getStartPos() );
						writer.write( "\t" + marked.getParser().getEndPos() );
						writer.write( "\t" + marked.getParser().getLength() );
						writer.write( "\t" + marked.getPercentIdentity());		
						writer.write( "\t" + marked.getShortestChainLength() + "\n");
						writer.flush();
					}	
				}
				
			}
		}	
		
		writer.close();							
	}	
	
}
