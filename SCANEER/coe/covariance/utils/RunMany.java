package covariance.utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.HashSet;
import java.util.List;

import scripts.PdbDownload;

import covariance.datacontainers.Alignment;
import covariance.datacontainers.AlignmentLine;
import covariance.datacontainers.PdbFileWrapper;
import covariance.parsers.BestAlignmentResultsFileLine;
import covariance.parsers.PfamParser;

public class RunMany
{
	// this will turn off filtering.
	// if set to 0 < filter < 100, sequences with percent identity 
	// above the filter percentage will be excluded from the alignemnts.
	public static final float FILTER_VAL = -1;
	
	public static void main(String[] args) throws Exception
	{
		PfamParser pFamParser = new PfamParser();
											
		List<BestAlignmentResultsFileLine> bestAlignmentList = 
			BestAlignmentResultsFileLine.getBestAlignmentResultsList(
					ConfigReader.getHomeDirectory() + File.separator + 
					"AlignmentScanResults.txt"
					);
		
		HashSet<String> pfamIdsToInclude = 
			BestAlignmentResultsFileLine.getPfamIdsAsSet(bestAlignmentList);
		
		while( ! pFamParser.isFinished()) 
		{
			Alignment a = pFamParser.getNextAlignment();
			
			if ( a != null)
			{
				System.out.println("Parsing " + a.getAligmentID());
				
				if ( pfamIdsToInclude.contains(a.getAligmentID()))
				{
					System.out.println("Analyzing " + a.getAligmentID() + " with " + 
							a.getNumSequencesInAlignment() + " sequences ");
					BestAlignmentResultsFileLine rFileLine = 
						BestAlignmentResultsFileLine.getOneById(bestAlignmentList, a.getAligmentID());
					
					File pdbFile =  PdbDownload.getOrDownloadFile(rFileLine.getPdbID());
					PdbFileWrapper pdbFileWrapper = new PdbFileWrapper(pdbFile);
					AlignmentLine aLine = a.getAnAlignmentLine(rFileLine.getAlignmentLineId());
		
					RefSequenceWrapper refSeqWrapper = 
							new RefSequenceWrapper(pdbFileWrapper,
									rFileLine, 
									a,
									aLine,
									FILTER_VAL);
						
					File outFile = new File( ConfigReader.getOutDataDir() + File.separator 
								+ a.getAligmentID() +  ".txt" );
						
					System.out.println("writing to " + outFile.getAbsolutePath() );
						
					BufferedWriter writer = new BufferedWriter( new FileWriter(outFile ));
					RunOne.writeFileHeader(writer);
					RunOne.runRefSequence(refSeqWrapper, writer);
				}
			}
		}
	}
}
