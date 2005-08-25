/*
 * History
 * Created on Sep 20, 2004 by kofke
 */
package etomica.nbr;

import etomica.atom.AtomPositionDefinition;
import etomica.atom.AtomSequencerFactory;
import etomica.atom.AtomType;
import etomica.atom.iterator.IteratorDirective;
import etomica.nbr.cell.IteratorFactoryCell;
import etomica.nbr.cell.NeighborCellManager;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.nbr.list.AtomSequencerNbr;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterNbr;
import etomica.phase.Phase;
import etomica.potential.Potential;
import etomica.potential.PotentialCalculation;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.species.Species;

/**
 * PotentialMaster that uses both neighbor-cell iteration and cell-list 
 * iteration.  This is needed by simulations that employ both Monte Carlo
 * and molecular dynamics integration steps, alternately as the simulation
 * proceeds.  See DCVGCMD simulation module for an example.
 * <br>
 */
public class PotentialMasterHybrid extends PotentialMaster {

	/**
	 * Invokes superclass constructor, specifying IteratorFactoryCell
     * for generating molecule iterators.  Sets default nCells of 10,
     * and position definition to null, causing cell assignment to be
     * based on atom type's position definition. 
	 */
	public PotentialMasterHybrid(Space space, double range) {
        this(space, null, range);
    }
    
    /**
     * Constructs class using given position definition for all atom cell assignments.
     * @param positionDefinition if null, specifies use of atom type's position definition
     */
   public PotentialMasterHybrid(Space space, AtomPositionDefinition positionDefinition, double range) {
        super(space,new IteratorFactoryCell());
        potentialMasterNbr = new PotentialMasterNbr(space, range, positionDefinition);
        potentialMasterCell = new PotentialMasterCell(space, range, positionDefinition);
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
     * Forward addPotentialBad to the PotentialMasterCell.  PotentialMasterNbr 
     * needs this too, but gets it on its own from NeighborManager.
     */
    public void addPotentialBad(Potential potential, AtomType[] atomTypes) {
        potentialMasterCell.addPotentialBad(potential, atomTypes);
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
    
    public NeighborListManager getNeighborManager() {
        return potentialMasterNbr.getNeighborManager();
    }

    private boolean useNbrLists;
    private final PotentialMasterNbr potentialMasterNbr;
    private final PotentialMasterCell potentialMasterCell;
}
