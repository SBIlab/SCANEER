package scripts;


import covariance.utils.ConfigReader;
import covariance.utils.ProcessWrapper;

import java.io.File;

import com.enterprisedt.net.ftp.FTPClient;
import com.enterprisedt.net.ftp.FTPTransferType;

public class PdbDownload
{
	public static final String URL_PREFIX="http://www.rcsb.org/pdb/cgi/export.cgi/";
	
	/**  Attempts to download pdbId from the pdb database.
	 *   Requires that ConfigReader can find a gzip binary 
	 */
	public static void downloadPdbFile( String pdbID) throws Exception
	{
		FTPClient ftpClient = new FTPClient("ftp.rcsb.org");
		ftpClient.login( "anonymous", "afodor@stanford.edu");
		ftpClient.setType(FTPTransferType.BINARY);
		
		if ( pdbID.length() != 4 )
			throw new Exception("Error!  Expecting pdbID for " + pdbID );
			
		ftpClient.get( getFilePath( pdbID) + ".Z", "/pub/pdb/data/structures/all/pdb/pdb" + pdbID + ".ent.Z" );
		ftpClient.quit();
				
		String[] args = new String[3];
		args[0] = ConfigReader.getGZipPath();
		args[1] = "-d";
		args[2] = getFilePath( pdbID) + ".Z";
		ProcessWrapper pw = new ProcessWrapper( args );
		
	}
	
	/**
	 *   If the file doesn't exist, attempt to get it off the internet.
	 *   Throws if this can't be done
	 *
	 */
	public static File getOrDownloadFile(String fourCharID) throws Exception
	{
		File pdbFile = PdbDownload.getFilePathFromString(fourCharID);
		
		if ( ! pdbFile.exists() || pdbFile.length() < 50)
		{
			System.out.println("Could not find " + pdbFile.getAbsolutePath() );
			System.out.println("Attempting download of " + fourCharID);
			PdbDownload.downloadPdbFile(fourCharID);
			if ( ! pdbFile.exists() || pdbFile.length() < 50)
				throw new Exception("Could not download " +  fourCharID);
		}
		else
		{
			System.out.println("Found " + pdbFile.getAbsolutePath());
		}
		
		return pdbFile;
	}
	
	public static File getFilePathFromString( String pdbId) throws Exception
	{
		return new File(getFilePath(pdbId));
	}
	
	public static String getFilePath( String pdbID ) throws Exception
	{
		return ConfigReader.getLocalPdbDirectory() + 
			   File.separator + pdbID ;
	}
}
