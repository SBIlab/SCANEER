package covariance.utils;

import java.io.*;
import java.util.*;

public class ConfigReader
{
	public static final String ENERGETICS_PROPERTIES_FILE="Energetics.properties";
	
	public static final String CLUSTAL_EXECUTABLE="CLUSTAL_EXECUTABLE";
	public static final String CLUSTAL_DIRECTORY="CLUSTAL_DIRECTORY";
	public static final String HOME_DIRECTORY="HOME_DIRECTORY";
	public static final String FULL_PFAM_PATH="FULL_PFAM_PATH";
	public static final String LOCAL_PDB_DIRECTORY="LOCAL_PDB_DIRECTORY";
	public static final String OUT_DATA_DIR="OUT_DATA_DIR";
	public static final String PY_MOL_OUT_DIR="PY_MOL_OUT_DIR";
	public static final String NATIVE_LOCKLESS_DIR="NATIVE_LOCKLESS_DIR";
	public static final String GZIP_FULL_PATH="GZIP_FULL_PATH";
	
	public static final String TRUE="TRUE";
	public static final String YES="YES";
	
	private static ConfigReader configReader = null;
	private static Properties props = new Properties();
	
	private static String getAProperty(String namedProperty ) throws Exception
	{
		Object obj = props.get( namedProperty );
		
		if ( obj == null ) 
			throw new Exception("Error!  Could not find " + namedProperty + " in " + ENERGETICS_PROPERTIES_FILE );
		
		return obj.toString();
	}
	
	private ConfigReader() throws Exception
	{
		Object o = new Object();
		
		InputStream in = o.getClass().getClassLoader().getSystemResourceAsStream( ENERGETICS_PROPERTIES_FILE );
		
		if ( in == null )
			throw new Exception("Error!  Could not find " + ENERGETICS_PROPERTIES_FILE + " anywhere on the current classpath");		
		
		BufferedReader reader = new BufferedReader( new InputStreamReader( in ));
		props = new Properties();
		
		String nextLine = reader.readLine();
		
		while ( nextLine != null ) 
		{
			nextLine= nextLine.trim();
			
			if ( nextLine.length() > 0 && ! nextLine.startsWith("!") && ! nextLine.startsWith("#") ) 
			{
				StringTokenizer sToken = new StringTokenizer( nextLine, "=" );
				
				if ( sToken.hasMoreTokens() ) 
				{
					String key = sToken.nextToken().trim();
					
					String value = "";
					
					if ( sToken.hasMoreTokens() ) 
					{
						value = sToken.nextToken().trim();
					}
					
					props.put( key, value );
					
				}	
			}	
			nextLine = reader.readLine();
		}	
	}
	
	private static synchronized ConfigReader getConfigReader() throws Exception
	{
		if ( configReader == null ) 
		{
			configReader = new ConfigReader();
		}
		
		return configReader;
	}
	
	public static String getHomeDirectory() throws Exception
	{
		return getConfigReader().getAProperty( HOME_DIRECTORY );
	}
	
	public static String getFullPfamPath() throws Exception
	{
		return getConfigReader().getAProperty( FULL_PFAM_PATH );
	}
	
	public static String getLocalPdbDirectory() throws Exception
	{
		return getConfigReader().getAProperty( LOCAL_PDB_DIRECTORY );
	}
	
	public static String getOutDataDir() throws Exception
	{
		return getConfigReader().getAProperty( OUT_DATA_DIR );
	}
	
	public static String getClustalExecutable() throws Exception
	{
		return getConfigReader().getAProperty( CLUSTAL_EXECUTABLE);
	}
	
	public static String getClustalDirectory() throws Exception
	{
		return getConfigReader().getAProperty( CLUSTAL_DIRECTORY);
	}
	
	public static String getPymolOutDir() throws Exception
	{
		return getConfigReader().getAProperty( PY_MOL_OUT_DIR);
	}
	
	public static String getNativeLocklessDir() throws Exception
	{
		return getConfigReader().getAProperty( NATIVE_LOCKLESS_DIR);
	}
	
	public static String getGZipPath() throws Exception
	{
		return getConfigReader().getAProperty(GZIP_FULL_PATH);
	}
}

