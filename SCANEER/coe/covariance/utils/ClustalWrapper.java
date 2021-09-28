package covariance.utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

import covariance.parsers.ClustalAlignment;

public class ClustalWrapper
{
	/** Not thread safe
	 */
	public static ClustalAlignment getClustalAlignment(String fasta) throws Exception
	{
		File clustalExecutable = new File( ConfigReader.getClustalExecutable());
		
		if ( ! clustalExecutable.isFile() || ! clustalExecutable.exists() ) 
			throw new Exception("Error!  Clustal file " + clustalExecutable.getAbsoluteFile() + " could not be found." +
								" Set " + ConfigReader.CLUSTAL_EXECUTABLE+ " in properties file");
		
		File fastaFile = writeTempFastaFile(fasta);
		
		String[] cmdArgs = new String[2];
		cmdArgs[0] = clustalExecutable.getAbsolutePath();
		cmdArgs[1] = fastaFile.getAbsolutePath();
		
		//System.out.println( cmdArgs[0] + " " + cmdArgs[1] );
		ProcessWrapper pw = new ProcessWrapper( cmdArgs );									
		
		//fastaFile.delete();	
		
		File alignmentFile = new File( ConfigReader.getClustalDirectory() + 
										File.separator + "temp.aln" );
		
		ClustalAlignment cA = new ClustalAlignment( alignmentFile );
		
		//alignmentFile.delete();
		
		return cA;
	}
	
	private static File writeTempFastaFile(String fasta) throws Exception
	{
		File fastaFile = new File( ConfigReader.getClustalDirectory()+ 
		 							File.separator + "temp" );	
		//System.out.println("Creating "+ fastaFile.getAbsolutePath() ); 						 							
		BufferedWriter writer = new BufferedWriter( new FileWriter( fastaFile ));
		writer.write(fasta);
		writer.flush();  writer.close();

		return fastaFile;								 							

	}
	
	/**  Just some test code
 	 */
	public static void main(String[] args) throws Exception
	{
		String fastaString = ">1\nELVIS\n>2\nELVIS\n";	
		System.out.println( getClustalAlignment(fastaString).getPairwiseIdentity(0,1));
	}
}
