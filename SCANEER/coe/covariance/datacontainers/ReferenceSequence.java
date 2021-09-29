package covariance.datacontainers;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;


public class ReferenceSequence
{
	private List referenceSequenceResidues = new ArrayList();
	
	public List getReferenceSequenceResidues()
	{
		return referenceSequenceResidues;
	}
	
	public void addReferenceSequenceResidue( ReferenceSequenceResidue refSequenceResidue ) 
	{
		this.referenceSequenceResidues.add( refSequenceResidue );
	}
	
	public ReferenceSequenceResidue getRefSeqByAlignmentPosition( int alignmentPosition) throws Exception
	{
		for ( Iterator i = referenceSequenceResidues.iterator();
				i.hasNext(); ) 
		{
			ReferenceSequenceResidue rSeqResidue = 
					(ReferenceSequenceResidue) i.next();
					
			if ( rSeqResidue.getAlignmentPosition() == alignmentPosition ) 
				return rSeqResidue;
		}
		
		return null;
	}
	
	public ReferenceSequenceResidue getRefSeqByPdbPosition( char pdbChainChar, int pdbPosition ) 
	{
		for ( Iterator i = referenceSequenceResidues.iterator();
				i.hasNext(); )
		{
			ReferenceSequenceResidue rSeqResidue = (ReferenceSequenceResidue)i.next();
			
			for ( Iterator i2= rSeqResidue.getLinkedPdbResidues().iterator();
					i2.hasNext(); ) 
			{
				PdbResidue pdbResidue = (PdbResidue) i2.next();
				
				if ( pdbResidue.getParentChain().getChainChar() == pdbChainChar && 
						pdbResidue.getPdbPosition() == pdbPosition ) 
						return rSeqResidue;
			}
			
		}
		
		return null;
	}
	
	public void calculateNeighbors(PdbFileWrapper pdbFileWrapper) throws Exception
	{
		for ( Iterator i = referenceSequenceResidues.iterator();
				i.hasNext(); ) 
		{
			ReferenceSequenceResidue rSeqResidue = ( ReferenceSequenceResidue ) i.next();
			
			if ( rSeqResidue.allMatches() ) 
			{
				for ( Iterator i2 = rSeqResidue.getLinkedPdbResidues().iterator();	
						i2.hasNext(); ) 
				{
					PdbResidue linkedResidue = ( PdbResidue )i2.next();
					linkedResidue.generateNeighborList( pdbFileWrapper.getPdbChains() );
				}
			}
		}	
	}
}
