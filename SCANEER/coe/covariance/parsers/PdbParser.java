package covariance.parsers;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.StringTokenizer;

import covariance.datacontainers.PdbAtom;
import covariance.datacontainers.PdbChain;
import covariance.datacontainers.PdbFileWrapper;
import covariance.datacontainers.PdbResidue;
import covariance.utils.SequenceUtils;


public class PdbParser
{
	private int numResidues = -1;
	
	public PdbParser(String filePathToParse, PdbFileWrapper pdbWrapper ) throws Exception
	{
		File fileToParse = new File(filePathToParse);
		
		if ( ! fileToParse.exists() ) 
			throw new Exception("Error!  Could not find " + fileToParse.getAbsolutePath() );
		
		List atomsList = readAtoms( fileToParse );
		populateWrapper(atomsList, pdbWrapper);
		pdbWrapper.setPdbId(getFourCharId(fileToParse).toLowerCase());
		pdbWrapper.setExperimentMethod(getExpData(fileToParse));
	}
	
	String getExpData( File fileToParse ) throws Exception
	{
		BufferedReader reader = new BufferedReader( new FileReader( fileToParse ));
		
		String nextLine = reader.readLine();
		
		while ( nextLine != null && nextLine.trim().length() > 0 )
		{
		
			if ( nextLine.startsWith("EXPDTA") ) 
			{
				StringTokenizer sToken = new StringTokenizer( nextLine );
				sToken.nextToken();
				
				return sToken.nextToken();
			}
			
			nextLine = reader.readLine();

		}
		
		return null;							
	}
	
	
	private String getFourCharId( File fileToParse ) throws Exception
	{
		BufferedReader reader = new BufferedReader( new FileReader( fileToParse ));
		
		String firstLine = reader.readLine();
		
		reader.close();
		
		return firstLine.substring(62,66);
	}
	
	private List readAtoms( File fileToParse ) throws Exception
	{
		List list = new ArrayList();
		
		BufferedReader reader = new BufferedReader( new FileReader( fileToParse));
		
		String nextLine = reader.readLine();
		
		while ( nextLine != null )
		{
			if ( nextLine.startsWith("ATOM"))
			{
				PdbAtom atom = new PdbAtom( nextLine );
				
				list.add( atom );
			}
			
			nextLine = reader.readLine();
		}		
		
		return list;
	}
	
	/**  Todo:  Implement more efficiently
	 */
	private HashSet getChainChars(List atomList)
	{
		HashSet chainChars = new HashSet();
		
		for (Iterator i = atomList.iterator();
			  i.hasNext(); ) 
		{
			PdbAtom atom = (PdbAtom) i.next();
			chainChars.add( new Character( atom.getChainId() ));
		}
		
		return chainChars;
	}
	
	/**  Todo:  Implement more efficiently
	 */
	private HashSet getResiduesInts(List atomList)
	{
		HashSet residueInts = new HashSet();
		
		for ( Iterator i = atomList.iterator();
			   i.hasNext(); ) 
		{
			PdbAtom pdbAtom = (PdbAtom) i.next();
			residueInts.add( new Integer( pdbAtom.getResidueSequenceNumber() ));
		}
		
		return residueInts;
	}
	
	private void populateWrapper( List atomList, PdbFileWrapper pdbWrapper ) throws Exception
	{
		Collection chainChars = getChainChars(atomList);
		
		for ( Iterator i = chainChars.iterator();
		       i.hasNext(); ) 
		{
			char chainChar = ((Character) i.next()).charValue();			
			PdbChain pdbChain = new PdbChain(chainChar, pdbWrapper);
			populateResidues(atomList, pdbChain);
			pdbWrapper.addChain(pdbChain);		
		}
	}
	
	private void populateResidues( List atomsList, PdbChain pdbChain ) throws Exception
	{
		Collection residueInts = getResiduesInts( atomsList);
		
		for ( Iterator i = residueInts.iterator();
				i.hasNext(); ) 
		{
			int residueInt = (( Integer ) i.next()).intValue();
			List subList = getAtomsInResidue(atomsList,  pdbChain.getChainChar(), residueInt);
			
			if ( subList.size() > 0 ) 
			{
				char residueChar = getResidueCharFromList(subList);
			
				PdbResidue pdbResidue = new PdbResidue(residueChar, residueInt, pdbChain);
				pdbChain.addResidue(pdbResidue);
				
				for ( Iterator i2 = subList.iterator();
					   i2.hasNext(); ) 
				{
					PdbAtom atom = (PdbAtom) i2.next();
					
					if ( pdbResidue.getAtom(atom.getAtomName()) != null ) 
						throw new Exception("Error!  Already have " + 
						 		atom.getResidueSequenceNumber() + " " +  atom.getAtomName() );					
						 		
					pdbResidue.addPdbAtom(atom);
				}
			}
		}
	}
	
	private char getResidueCharFromList( List subList ) throws Exception
	{
		String residueName = ((PdbAtom) subList.get(0)).getResidueName();
		
		for ( Iterator i = subList.iterator();
			   i.hasNext(); ) 
		{
			PdbAtom atom = (PdbAtom) i.next();
			
			if ( ! atom.getResidueName().equals(residueName) )
				throw new Exception("Error in pdb file.  Unequal residue names " + 
				  atom.getResidueName() + " vs " + residueName );
		}
		
		return SequenceUtils.threeToOne(residueName).charAt(0);
	}
	
	/**  Todo:  This is probably too slow
	 */
	private List getAtomsInResidue( List atomsList, 
									 char chainChar,
									 int residueInt ) throws Exception
	{
		List list = new ArrayList();
		
		for ( Iterator i = atomsList.iterator();
			   i.hasNext(); ) 
		{
			PdbAtom atom = (PdbAtom) i.next();
			
			if ( atom.getResidueSequenceNumber() == residueInt 
				  && atom.getChainId() == chainChar ) 
				list.add( atom );
		}
		
		return list;		
	}
}
