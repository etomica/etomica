/*
 * History
 * Created on Sep 20, 2004 by kofke
 */
package etomica.nbr;

import etomica.atom.AtomPositionDefinition;
import etomica.atom.AtomType;
import etomica.atom.iterator.IteratorDirective;
import etomica.nbr.cell.NeighborCellManager;
import etomica.nbr.cell.PhaseAgentSourceCellManager;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.phase.Phase;
import etomica.phase.PhaseAgentManager;
import etomica.potential.Potential;
import etomica.potential.PotentialCalculation;
import etomica.potential.PotentialGroup;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.species.Species;

/**
 * PotentialMaster that uses both neighbor-cell iteration and cell-list 
 * iteration.  This is needed by simulations that employ both Monte Carlo
 * and molecular dynamics integration steps, alternately as the simulation
 * proceeds.  See DCVGCMD simulation module for an example.
 * <br>
 */
public class PotentialMasterHybrid extends PotentialMasterNbr {

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
        this(space, range, new PhaseAgentSourceCellManager(positionDefinition));
    }
    
    private PotentialMasterHybrid(Space space, double range, PhaseAgentSourceCellManager phaseAgentSource) {
        this(space, range, phaseAgentSource, new PhaseAgentManager(phaseAgentSource,null));
    }
    
    private PotentialMasterHybrid(Space space, double range, PhaseAgentSourceCellManager phaseAgentSource,
            PhaseAgentManager agentManager) {
        super(space, phaseAgentSource, agentManager);
        potentialMasterList = new PotentialMasterList(space, range, phaseAgentSource, agentManager);
        potentialMasterCell = new PotentialMasterCell(space, range, phaseAgentSource, agentManager);
	}
    
    public PotentialGroup makePotentialGroup(int nBody) {
        return new PotentialGroupHybrid(nBody,space);
    }
    
    public void setSimulation(Simulation simulation) {
        super.setSimulation(simulation);
        potentialMasterList.setSimulation(simulation);
        potentialMasterCell.setSimulation(simulation);
    }
    
    public PotentialMasterList getPotentialMasterList() {
        return potentialMasterList;
    }
    
    public PotentialMasterCell getPotentialMasterCell() {
        return potentialMasterCell;
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
        if (useNbrLists) potentialMasterList.calculate(phase,id,pc);
        else potentialMasterCell.calculate(phase,id,pc);
    }
    
    public double getCellRange() {
        return potentialMasterCell.getRange();
    }
    
    public void setCellRange(int newRange) {
        potentialMasterList.setCellRange(newRange);
        potentialMasterCell.setCellRange(newRange);
    }
    
    public double getRange() {
        return potentialMasterCell.getRange();
    }

    public void setRange(double newRange) {
        potentialMasterList.setRange(newRange);
        potentialMasterCell.setRange(newRange);
    }

    public NeighborCellManager getNbrCellManager(Phase phase) {
        return potentialMasterList.getNbrCellManager(phase);
    }
    
    public void setUseNbrLists(boolean flag) {
        useNbrLists = flag;
    }
    
    public void addPotential(Potential potential, Species[] species) {
        potentialMasterList.addPotential(potential, species);
        potentialMasterCell.addPotential(potential, species);
        if (potential instanceof PotentialGroup) {
            // potential masters will attempt to set themselves as the group's 
            // PotentialMaster, but it will resist because it only has eyes for 
            // us.
            ((PotentialGroup)potential).setPotentialMaster(this);
        }
    }
    
    protected void addRangedPotentialForTypes(Potential potential, AtomType[] atomTypes) {
    }
    
    public void potentialAddedNotify(Potential subPotential, PotentialGroup pGroup) {
        potentialMasterList.potentialAddedNotify(subPotential, pGroup);
        potentialMasterCell.potentialAddedNotify(subPotential, pGroup);
    }
    
    public NeighborListManager getNeighborManager() {
        return potentialMasterList.getNeighborManager();
    }

    public void removePotential(Potential potential) {
        potentialMasterList.removePotential(potential);
        potentialMasterCell.removePotential(potential);
    }

    private boolean useNbrLists;
    private final PotentialMasterList potentialMasterList;
    private final PotentialMasterCell potentialMasterCell;
}
