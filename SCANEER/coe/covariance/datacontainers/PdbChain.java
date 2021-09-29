package covariance.datacontainers;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;

public class PdbChain
{
	private char chainChar;
	private HashSet pdbResidues = new HashSet();
	private PdbFileWrapper pdbFile;
	private int highestPdbResiduePosition = -1;
	
	public PdbChain( char chainId, PdbFileWrapper pdbFile ) 
	{
		this.chainChar = chainId;
		this.pdbFile = pdbFile;	
	}
	
	public void addResidue(PdbResidue pdbResidue)
	{
		this.pdbResidues.add( pdbResidue );	
		
		if( highestPdbResiduePosition < pdbResidue.getPdbPosition() ) 
			highestPdbResiduePosition = pdbResidue.getPdbPosition();
	}
	
	public String getSequence()
	{
		StringBuffer buff = new StringBuffer();
		
		List aList = new ArrayList( pdbResidues);
		
		Collections.sort( aList );
		
		for ( int x=0; x< aList.size(); x++ ) 
		{
			PdbResidue reside = (PdbResidue) aList.get(x);
			buff.append( reside.getPdbChar());
		}
		
		return buff.toString();
	}
	
	public char getChainChar()
	{
		return chainChar;
	}

	public PdbFileWrapper getPdbFile()
	{
		return pdbFile;
	}

	public HashSet getPdbResidues()
	{
		return pdbResidues;
	}
	
	
	/**  Todo:  Make more efficient
	 */
	public PdbResidue getPdbResidueByPdbPosition( int pdbPosition ) 
	{
		for ( Iterator i = pdbResidues.iterator();
			   i.hasNext(); ) 
		{
			PdbResidue pdbResidue = (PdbResidue) i.next();
			
			if ( pdbResidue.getPdbPosition() == pdbPosition ) 
				return pdbResidue;
		}
		
		return null;
	}

	public int getHighestPdbResiduePosition()
	{
		return highestPdbResiduePosition;
	}

}
