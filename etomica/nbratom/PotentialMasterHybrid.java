/*
 * History
 * Created on Sep 20, 2004 by kofke
 */
package etomica.nbratom;

import etomica.IteratorDirective;
import etomica.Phase;
import etomica.Potential;
import etomica.PotentialMaster;
import etomica.Space;
import etomica.Species;
import etomica.atom.AtomPositionDefinition;
import etomica.atom.AtomSequencerFactory;
import etomica.nbratom.cell.IteratorFactoryCell;
import etomica.nbratom.cell.NeighborCellManager;
import etomica.nbratom.cell.PotentialCalculationAgents;
import etomica.nbratom.cell.PotentialMasterCell;
import etomica.potential.PotentialCalculation;

/**
 * PotentialMaster used to implement neighbor listing.  Instance of this
 * class is given as an argument to the Simulation constructor.
 * Criteria specifying whether two atoms are neighbors for a particular potential
 * are specified in the setSpecies method of this class.
 * <br>
 */
public class PotentialMasterHybrid extends PotentialMaster {

	/**
	 * Invokes superclass constructor, specifying IteratorFactoryCell
     * for generating molecule iterators.  Sets default nCells of 10,
     * and position definition to null, causing cell assignment to be
     * based on atom type's position definition. 
	 */
	public PotentialMasterHybrid(Space space) {
        this(space, null);
    }
    
    /**
     * Constructs class using given position definition for all atom cell assignments.
     * @param positionDefinition if null, specifies use of atom type's position definition
     */
   public PotentialMasterHybrid(Space space, AtomPositionDefinition positionDefinition) {
        super(space,new IteratorFactoryCell());
        potentialMasterNbr = new PotentialMasterNbr(space, positionDefinition);
        potentialMasterCell = new PotentialMasterCell(space, positionDefinition);
	}
    
    /**
     * Performs cell-assignment potentialCalculation.  Assigns all molecules
     * to their cells, and invokes superclass method causing setup to be
     * performed iterating using species/potential hierarchy.
     */
    //TODO move this method to PotentialMaster superclass
    public void calculate(Phase phase, PotentialCalculationAgents pc) {
        super.calculate(phase,new IteratorDirective(),pc);
    }
    
    /**
     * Overrides superclass method to enable direct neighbor-list iteration
     * instead of iteration via species/potential hierarchy. If no target atoms are
     * specified in directive, neighborlist iteration is begun with
     * speciesMaster of phase, and repeated recursively down species hierarchy;
     * if one atom is specified, neighborlist iteration is performed on it and
     * down species hierarchy from it; if two or more atoms are specified,
     * superclass method is invoked.
     */
    public void calculate(Phase phase, IteratorDirective id, PotentialCalculation pc) {
		if(!enabled) return;
        if (useNbrLists) potentialMasterNbr.calculate(phase,id,pc);
        else potentialMasterCell.calculate(phase,id,pc);
	}//end calculate
	
    public int getNCells() {
        return potentialMasterCell.getNCells();
    }
    public void setNCells(int cells) {
        potentialMasterNbr.setNCells(cells);
        potentialMasterCell.setNCells(cells);
    }
    
    public NeighborCellManager getNbrCellManager(Phase phase) {
        NeighborCellManager manager = (NeighborCellManager)phase.getCellManager();
        if (manager == null) {
            manager = new NeighborCellManager(phase,getNCells(),getAtomPositionDefinition());
            phase.setCellManager(manager);
        }
        return manager;
    }

    public AtomSequencerFactory sequencerFactory() {return AtomSequencerNbr.FACTORY;}
    
    public AtomPositionDefinition getAtomPositionDefinition() {
        return potentialMasterNbr.getAtomPositionDefinition();
    }

    public void setUseNbrLists(boolean flag) {
        useNbrLists = flag;
    }
    
    /* (non-Javadoc)
     * @see etomica.PotentialMaster#setSpecies(etomica.Potential, etomica.Species[])
     */
    public void setSpecies(Potential potential, Species[] species) {
        potentialMasterNbr.setSpecies(potential, species);
        potentialMasterCell.setSpecies(potential, species);
        super.setSpecies(potential, species);
    }
    
    public NeighborManager getNeighborManager() {
        return potentialMasterNbr.getNeighborManager();
    }

    private boolean useNbrLists;
    private final PotentialMasterNbr potentialMasterNbr;
    private final PotentialMasterCell potentialMasterCell;
}
