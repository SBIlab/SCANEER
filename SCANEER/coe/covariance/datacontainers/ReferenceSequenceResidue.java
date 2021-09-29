package covariance.datacontainers;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class ReferenceSequenceResidue
{
	private char referenceChar;
	private int alignmentPosition;
	private List linkedPdbResidues = new ArrayList();
	
	public boolean allMatches()
	{
		if ( referenceChar == Alignment.GAP_CHAR ) 
			return false;
			
		for ( Iterator i = linkedPdbResidues.iterator();
				i.hasNext(); ) 
		{
			PdbResidue pdbResidue = (PdbResidue) i.next();
			if ( pdbResidue.getPdbChar() != referenceChar  || pdbResidue.getCbAtom() == null  )
				return false;
		}		
		
		return true;
	}	
	
	public int getAlignmentPosition()
	{
		return alignmentPosition;
	}

	public List getLinkedPdbResidues()
	{
		return linkedPdbResidues;
	}

	public char getReferenceChar()
	{
		return referenceChar;
	}
	
	public PdbResidue getLinkedPdbResidue( char chainChar)
	{
		for ( Iterator i = this.linkedPdbResidues.iterator();
				i.hasNext(); )	
		{
			PdbResidue pdbRes = (PdbResidue) i.next();
			
			if ( pdbRes.getParentChain().getChainChar() == chainChar ) 
				return pdbRes;
		}
		
		return null;
	} 
	
	
	public double getMinCbDistance( ReferenceSequenceResidue otherResidue ) throws Exception
	{
		double min = Double.MAX_VALUE;
		
		for ( Iterator i = this.linkedPdbResidues.iterator();
				i.hasNext(); ) 
		{
			PdbResidue iResidue = (PdbResidue) i.next();
			PdbAtom iAtom = iResidue.getCbAtom();
			
			for ( Iterator j = otherResidue.linkedPdbResidues.iterator();
						j.hasNext(); ) 
			{
				PdbResidue jResidue = (PdbResidue) j.next();
				PdbAtom jAtom = jResidue.getCbAtom();	
				
				double score = iAtom.getDistance( jAtom);
						
				if ( score < min ) 
					min = score;
			}
		}
		
		if ( min == Double.MAX_VALUE ) 
			throw new Exception("Comparison not valid");
		
		return min;
	}
	
	public double getAvgNumAtomAtomContacts( List pdbChains ) throws Exception
	{
		double sum = 0;
		int n =0;
		
		for ( Iterator i = this.getLinkedPdbResidues().iterator(); 
				i.hasNext(); ) 
		{
			PdbResidue pdbResidue = (PdbResidue) i.next();
			sum += pdbResidue.getTotalNumberAtomAtomContactsFromChains( pdbChains);
			n++;
		}
		
		return sum / n;	
	}
	
	public double getAvgNumContacts( ReferenceSequenceResidue otherResidue ) throws Exception
	{
		double sum = 0;
		int n = 0;
		
		for ( Iterator i = this.getLinkedPdbResidues().iterator();
				i.hasNext(); ) 
		{
			PdbResidue pdbResidue = (PdbResidue ) i.next();
			sum += pdbResidue.getNumNeighbors();
			n++;
		}
		
		for ( Iterator i = otherResidue.getLinkedPdbResidues().iterator();
				i.hasNext(); ) 
		{
			PdbResidue pdbResidue = (PdbResidue ) i.next();
			sum += pdbResidue.getNumNeighbors();
			n++;
		}
		
		return sum/n;		
	
	}
	
	public double getAvgNumContacts() throws Exception
	{
		double sum = 0;
		int n = 0;
		
		for ( Iterator i = this.getLinkedPdbResidues().iterator();
				i.hasNext(); ) 
		{
			PdbResidue pdbResidue = (PdbResidue ) i.next();
			sum += pdbResidue.getNumNeighbors();
			n++;
		}
		
		return sum/n;		
	
	}
	
	public double getAvgNumHydrophobicContacts() throws Exception
	{
		double sum = 0;
		int n = 0;
		
		for ( Iterator i = this.getLinkedPdbResidues().iterator();
				i.hasNext(); ) 
		{
			PdbResidue pdbResidue = (PdbResidue ) i.next();
			sum += pdbResidue.getNumHydrophicNeighbors();
			n++;
		}
	
		return sum/n;
	}
	
	public double getAvgNumHydrophobicContacts( ReferenceSequenceResidue otherResidue ) throws Exception
	{
		double sum = 0;
		int n = 0;
		
		for ( Iterator i = this.getLinkedPdbResidues().iterator();
				i.hasNext(); ) 
		{
			PdbResidue pdbResidue = (PdbResidue ) i.next();
			sum += pdbResidue.getNumHydrophicNeighbors();
			n++;
		}
		
		for ( Iterator i = otherResidue.getLinkedPdbResidues().iterator();
				i.hasNext(); ) 
		{
			PdbResidue pdbResidue = (PdbResidue ) i.next();
			sum += pdbResidue.getNumHydrophicNeighbors();
			n++;
		}
		
		return sum/n;		
	
	}
	
	public double getAvgTemperature()
	{
		double sum = 0;
		int n = 0;
		
		List atomList = new ArrayList();
		
		for ( Iterator i = this.getLinkedPdbResidues().iterator();
				i.hasNext(); ) 
		{
			PdbResidue pdbResidue = (PdbResidue ) i.next();
			atomList.addAll( pdbResidue.getPdbAtoms() );
		}
		
		for ( Iterator i = atomList.iterator();
				i.hasNext(); ) 
		{
			PdbAtom atom = (PdbAtom) i.next();
			sum += atom.getTemperature();	
			n++;
		}
		
		return sum/n;		
		
	}
	
	public double getAvgTemperature( ReferenceSequenceResidue otherResidue ) 
	{
		double sum = 0;
		int n = 0;
		
		List atomList = new ArrayList();
		
		for ( Iterator i = this.getLinkedPdbResidues().iterator();
				i.hasNext(); ) 
		{
			PdbResidue pdbResidue = (PdbResidue ) i.next();
			atomList.addAll( pdbResidue.getPdbAtoms() );
		}
		
		for ( Iterator i = otherResidue.getLinkedPdbResidues().iterator();
				i.hasNext(); ) 
		{
			PdbResidue pdbResidue = (PdbResidue ) i.next();
			atomList.addAll( pdbResidue.getPdbAtoms() );
		}
		
		for ( Iterator i = atomList.iterator();
				i.hasNext(); ) 
		{
			PdbAtom atom = (PdbAtom) i.next();
			sum += atom.getTemperature();	
			n++;
		}
		
		return sum/n;		
	}
	
	public ReferenceSequenceResidue( char referenceChar,
									  int alignmentPosition ) 
	{
		this.referenceChar = referenceChar;
		this.alignmentPosition = alignmentPosition;
	}									  	

	public void addLinkedPdbResidue( PdbResidue pdbResidue ) 
	{
		this.linkedPdbResidues.add(pdbResidue);
	}
}
