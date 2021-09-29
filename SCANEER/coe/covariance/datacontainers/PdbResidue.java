package covariance.datacontainers;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import covariance.utils.SequenceUtils;

public class PdbResidue implements Comparable
{
	private char pdbChar;
	private int pdbPosition;
	private PdbChain parentChain;
	private List pdbAtoms = new ArrayList();
	private int numNeighbors = -1;
	private int numHydrophobicNeighbors = -1;
	private PdbAtom cbAtom;
	private boolean isHydrophobic;
	
	public static final float NEIGHBOR_DISTANCE_CUTOFF=10;
	public static final char BLANK_CHAR='-';
	public static final int BLANK_POSITION=-99;
	
	public PdbChain getParentChain()
	{
		return parentChain;
	}
	
	
	
	public boolean isHydrophobic()
	{
		return isHydrophobic;
	}

	public void generateNeighborList(List pdbChains)
	{
		this.numNeighbors = 0;
		this.numHydrophobicNeighbors = 0;
		
		for ( Iterator i = pdbChains.iterator();
				i.hasNext(); )
		{
			PdbChain chain = (PdbChain) i.next();
			
			for ( Iterator i2 = chain.getPdbResidues().iterator();
					i2.hasNext(); ) 
			{
				PdbResidue otherResidue = (PdbResidue) i2.next();
				if ( otherResidue.cbAtom != null &&
						this != otherResidue &&
						this.cbAtom.getDistance( otherResidue.cbAtom ) < NEIGHBOR_DISTANCE_CUTOFF ) 
				{
					this.numNeighbors++;
					
					if ( otherResidue.isHydrophobic ) 
						this.numHydrophobicNeighbors++;			
				}
			}
		}	
	}
	
	/** call generateNeighborList before calling this method.   
	 */
	public int getNumNeighbors() throws Exception
	{
		if ( this.numNeighbors == -1 ) 
			throw new Exception("Error!  call generateNeighborList() before calling this method" );
		
		return this.numNeighbors;
	}
	
	/** call generateNeighborList before calling this method. 
	 */	
	public int getNumHydrophicNeighbors() throws Exception
	{
		if ( this.numHydrophobicNeighbors== -1 ) 
			throw new Exception("Error!  call generateNeighborList() before calling this method" );
		
		return this.numHydrophobicNeighbors;
		
	}

	public char getPdbChar()
	{
		return pdbChar;
	}

	public int getPdbPosition()
	{
		return pdbPosition;
	}
	
	public void addPdbAtom( PdbAtom pdbAtom ) 
	{
		if ( pdbChar == 'G' ) 
		{ 
			if (pdbAtom.getAtomName().equals("CA") ) 
				this.cbAtom = pdbAtom;
		}
		else
		{
			if ( pdbAtom.getAtomName().equals("CB") ) 
				this.cbAtom = pdbAtom;
		}
		
		this.pdbAtoms.add( pdbAtom );	
	}
	
	public PdbResidue( char pdbChar, int pdbPosition, PdbChain parentChain ) 
	{
		this.pdbChar = pdbChar;
		this.pdbPosition = pdbPosition;
		this.parentChain = parentChain;	
		
		if (SequenceUtils.isHydrophicChar( this.pdbChar) )
		{
			this.isHydrophobic = true;	
		}
		else
		{
			this.isHydrophobic = false;
		}
	}
	
	public int getTotalNumberAtomAtomContactsFromChains( List pdbChains ) throws Exception
	{
		int sum = 0;
		
		for ( Iterator i = pdbChains.iterator();
				 i.hasNext(); )
		{
			PdbChain chain = (PdbChain) i.next();
			sum += getTotalNumberAtomAtomContacts( chain.getPdbResidues() );
		}
		
		return sum;		 
	}
	
	private int getTotalNumberAtomAtomContacts( Collection pdbResidues ) throws Exception
	{
		int sum = 0;
		
		for ( Iterator i = pdbResidues.iterator();
				i.hasNext(); ) 
		{
			PdbResidue pdbResidue = (PdbResidue) i.next();
			
			if ( this != pdbResidue )
			{
				sum+= getNumAtomAtomContacts(pdbResidue);
			}
		}
		
		return sum;	
	}
	
	/**  Returns the number of atoms within 5A of each other between this residue and otherResidue
	 */
	private int getNumAtomAtomContacts( PdbResidue otherResidue ) throws Exception
	{
		int returnInt =0;
		
		for ( Iterator i1 = this.pdbAtoms.iterator();
			   i1.hasNext(); ) 
		{
			PdbAtom thisAtom = (PdbAtom) i1.next();
								
			for (Iterator i2= otherResidue.pdbAtoms.iterator();
			   		 i2.hasNext(); ) 
			{
				PdbAtom otherAtom = (PdbAtom) i2.next();
				
				if ( thisAtom.getDistance( otherAtom) <= 5.0f ) 		   	
					returnInt++;	
			}
		}
		
		return returnInt;		
	}
	public float getAverageTemperature() throws Exception
	{
		if ( this.pdbAtoms.size() == 0 ) 
			throw new Exception("Error!  No atoms in residue ");
		
		float sum =0;
		
		for ( Iterator i = this.pdbAtoms.iterator();
				i.hasNext(); )
		{
			PdbAtom atom = (PdbAtom) i.next();
			sum += atom.getTemperature();		
		}
		
		return sum / this.pdbAtoms.size();
	}
	
	/**  todo:  Make more efficient
	 * 
	 *   returns null if the atom can't be found   
	 */
	public PdbAtom getAtom( String atomName ) 
	{
		for ( Iterator i = pdbAtoms.iterator();
			   i.hasNext();
			 )
		{
			PdbAtom atom = (PdbAtom) i.next();
			if ( atom.getAtomName().equals(atomName) ) 
				return atom;			
		}
	
		return null;	
	}

	public PdbAtom getCbAtom()
	{
		return this.cbAtom;
	}
	
	
	public List getPdbAtoms()
	{
		return pdbAtoms;
	}

	public int compareTo(Object o)
	{
		PdbResidue other = (PdbResidue ) o;
	
		return Float.compare(this.pdbPosition, other.pdbPosition);
	}

}
