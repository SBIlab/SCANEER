package covariance.utils;

import java.util.ArrayList;
import java.util.List;

import covariance.datacontainers.Alignment;
import covariance.datacontainers.AlignmentLine;
import covariance.datacontainers.PdbAtom;
import covariance.datacontainers.PdbChain;
import covariance.datacontainers.PdbFileWrapper;
import covariance.datacontainers.PdbResidue;
import covariance.datacontainers.ReferenceSequence;
import covariance.datacontainers.ReferenceSequenceResidue;
import covariance.parsers.BestAlignmentResultsFileLine;
import covariance.parsers.ClustalAlignment;

public class RefSequenceWrapper
{
	private final PdbFileWrapper pdbFileWrapper;
	private final PdbChain pdbChain;
	private final Alignment unFilteredAlignment;
	private final Alignment filteredAlignment;
	private final int pdbStart;
	private final int pdbEnd;
	private final char chainChar;
	private final AlignmentLine alignmentLine;
	private final BestAlignmentResultsFileLine rFileLine;
	private final ClustalAlignment cAlignment;
	private final ReferenceSequence referenceSequence;
	private final float filter;
	
	public RefSequenceWrapper( PdbFileWrapper pdbFileWrapper,
								BestAlignmentResultsFileLine rFileLine, 
								Alignment a,
								AlignmentLine alignmentLine,
								float filter) throws Exception
	{
		this(pdbFileWrapper, 
			rFileLine.getPdbChain(),
			rFileLine.getPdbStartResidue(),
			rFileLine.getPdbEndResidue(),
			a,
			alignmentLine,
			filter
			);
	}
	
	/**  set the filter to <= 0 to disable filtering
	 */
	public RefSequenceWrapper( PdbFileWrapper pdbFileWrapper,
								char chainChar,
								int pdbStart,
								int pdbEnd,
								Alignment alignment,
								AlignmentLine aLine,
								float filter ) throws Exception
	{
		this.pdbFileWrapper = pdbFileWrapper;
		this.pdbChain = pdbFileWrapper.getChain(chainChar);
		this.unFilteredAlignment= alignment;
		this.filter = filter;
		
		if ( filter >= 100 ) 
			throw new Exception("Error!  Set filter to <=0 to disable filtering.  Otherwise, filter should " + 
					" be 0 < filter < 100 ");
		
		if ( filter <= 0 ) 
		{
			this.filteredAlignment = alignment;
		}
		else
		{
			System.out.println("Removing all sequences > " + filter + " redundant starting with " + 
							alignment.getNumSequencesInAlignment() + " sequences " );
			System.out.println("This may be slow");
			this.filteredAlignment = alignment.getFilteredAlignment(filter);
			System.out.println("Finished with " + this.filteredAlignment.getNumSequencesInAlignment() );	
		}
	
		List list = new ArrayList();
		list.add( "#=GF DR   PDB; " + pdbFileWrapper.getFourCharId() + " " + chainChar
										+ "; " + pdbStart +"; " +pdbEnd+";");
		
		unFilteredAlignment.setPdbIds(list);
		filteredAlignment.setPdbIds(list);
		
		this.pdbStart = pdbStart;
		this.pdbEnd = pdbEnd;
		this.chainChar = chainChar;
		this.alignmentLine = aLine;
		
		this.rFileLine = new BestAlignmentResultsFileLine( unFilteredAlignment, alignmentLine, 
				pdbFileWrapper.getFourCharId(),
				chainChar, pdbStart, pdbEnd);
				
		this.cAlignment = ClustalWrapper.getClustalAlignment( 
		ReferenceSequenceUtils.buildFasta(unFilteredAlignment, alignmentLine, rFileLine, pdbFileWrapper));
		
		this.referenceSequence= ReferenceSequenceUtils.getReferenceSequence(unFilteredAlignment, alignmentLine, 
								cAlignment, rFileLine, pdbFileWrapper);
								
		referenceSequence.calculateNeighbors(pdbFileWrapper);
		
	}
	
	public void dumpRefSeqToConsole()
	{
		List residues = this.referenceSequence.getReferenceSequenceResidues();
		System.out.println( "alignmentPos\tpdbPos\tpdbChar\n");
		
		for ( int x=0; x < residues.size(); x++ ) 
		{
			ReferenceSequenceResidue refSeqResidue = (ReferenceSequenceResidue) residues.get(x);	
			PdbResidue pdbResidue = refSeqResidue.getLinkedPdbResidue(this.chainChar);
			
			System.out.println( refSeqResidue.getAlignmentPosition() + "\t" + 
										pdbResidue.getPdbPosition() + "\t" + 
											pdbResidue.getPdbChar() );
		}
	}
		
	
	public PdbAtom getCbAtom( List list, int listPosition ) throws Exception
	{
		ReferenceSequenceResidue refSeqResidue = 
					 (ReferenceSequenceResidue)	list.get(listPosition);
					 
		if ( refSeqResidue == null ) 
			return null;
			
		PdbResidue pdbResidue = refSeqResidue.getLinkedPdbResidue(this.chainChar);
		
		if ( pdbResidue == null ) 
			return null;
			
		return pdbResidue.getCbAtom();
		
	}
	
	public float getAverageCbDistance() throws Exception
	{
		float sum = 0;
		int n =0;
		
		List list = this.referenceSequence.getReferenceSequenceResidues();
		
		for ( int x=0; x < list.size(); x++ ) 
		{
			PdbAtom xAtom = getCbAtom(list, x);
			
			if ( xAtom != null  )
			{
				for ( int y=0; y < list.size(); y++ ) 
				{
					PdbAtom yAtom = getCbAtom(list, y);
					
					if ( yAtom != null ) 
					{
						sum+= xAtom.getDistance(yAtom);
						n++;
					}		
				}
			}	 
		}		
		
		
		return sum/n;
	}
	
	
	public AlignmentLine getAlignmentLine()
	{
		return alignmentLine;
	}

	public ClustalAlignment getCAlignment()
	{
		return cAlignment;
	}

	public char getChainChar()
	{
		return chainChar;
	}

	public PdbChain getPdbChain()
	{
		return pdbChain;
	}

	public int getPdbEnd()
	{
		return pdbEnd;
	}

	public PdbFileWrapper getPdbFileWrapper()
	{
		return pdbFileWrapper;
	}

	public int getPdbStart()
	{
		return pdbStart;
	}

	public ReferenceSequence getReferenceSequence()
	{
		return referenceSequence;
	}

	public BestAlignmentResultsFileLine getRFileLine()
	{
		return rFileLine;
	}
	
	public float getFilter()
	{
		return filter;
	}

	public Alignment getFilteredAlignment()
	{
		return filteredAlignment;
	}

	public Alignment getUnFilteredAlignment()
	{
		return unFilteredAlignment;
	}

}
