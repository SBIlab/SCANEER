package covariance.datacontainers;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import covariance.parsers.PdbParser;
import covariance.utils.ConfigReader;

public class PdbFileWrapper
{
	private String pdbId;
	private List pdbChains = new ArrayList();
	
	private String experimentMethod;
	
	
	
	public PdbFileWrapper( File pdbFile ) throws Exception
	{
		new PdbParser(pdbFile.getAbsolutePath() , this);
	}
	
	public PdbFileWrapper( String fourCharId ) throws Exception
	{
		this(new File( ConfigReader.getLocalPdbDirectory() + 
							File.separator + fourCharId) );
	}
	
	public List getPdbChains()
	{
		return pdbChains;
	}

	public String getFourCharId()
	{
		return pdbId;
	}

	public void setPdbId(String pdbId)
	{
		this.pdbId = pdbId;
	}
	
	public void addChain(PdbChain pdbChain ) 
	{
		pdbChains.add(pdbChain);
	}
	
	public PdbChain getChain( Character chainChar ) 
	{
		return getChain( chainChar.charValue());
	}
	
	public PdbChain getChain( char chainChar ) 
	{
		for ( Iterator i = this.pdbChains.iterator();
					i.hasNext(); ) 
		{
			PdbChain pdbChain = ( PdbChain ) i.next();
			
			if ( pdbChain.getChainChar() == chainChar ) 
				return pdbChain;
		}
		
		return null;
	}
	
	public int getLongestLength() 
	{
		int longestLength = -1;
		
		for ( int x=0; x < this.pdbChains.size(); x++ ) 
		{
			PdbChain chain= (PdbChain) this.pdbChains.get(x);
			
			longestLength = Math.max( longestLength, chain.getPdbResidues().size() );	
		}
		
		return longestLength;
	}

	public String getExperimentMethod()
	{
		return experimentMethod;
	}

	public void setExperimentMethod(String experimentMethod)
	{
		this.experimentMethod = experimentMethod;
	}

}
