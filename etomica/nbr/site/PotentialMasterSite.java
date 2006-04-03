/*
 * History
 * Created on Sep 20, 2004 by kofke
 */
package etomica.nbr.site;

import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.atom.iterator.AtomsetIteratorMolecule;
import etomica.atom.iterator.AtomsetIteratorSinglet;
import etomica.atom.iterator.IteratorDirective;
import etomica.nbr.NeighborCriterion;
import etomica.nbr.PotentialMasterNbr;
import etomica.nbr.PotentialCalculationUpdateTypeList.PotentialAtomTypeWrapper;
import etomica.nbr.cell.NeighborCellManager;
import etomica.phase.Phase;
import etomica.phase.PhaseAgentManager;
import etomica.phase.PhaseCellManager;
import etomica.phase.PhaseAgentManager.PhaseAgentSource;
import etomica.potential.Potential;
import etomica.potential.Potential2;
import etomica.potential.PotentialCalculation;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.util.Arrays;

public class PotentialMasterSite extends PotentialMasterNbr {

	/**
	 * Invokes superclass constructor, specifying IteratorFactoryCell
     * for generating molecule iterators.  Sets default nCells of 10 and
     * position definition to null, so that atom type's definition is used
     * to assign cells. 
	 */
	public PotentialMasterSite(Space space, int nCells) {
        this(space, new PhaseAgentSiteManager(nCells));
    }
    
    public PotentialMasterSite(Space space, PhaseAgentSource phaseAgentSource) {
        this(space, phaseAgentSource, new PhaseAgentManager(phaseAgentSource,null));
    }
    
    public PotentialMasterSite(Space space, PhaseAgentSource phaseAgentSource, PhaseAgentManager agentManager) {
        this(space, phaseAgentSource, agentManager, new Api1ASite(space.D(),agentManager));
    }
    
    protected PotentialMasterSite(Space space, PhaseAgentSource phaseAgentSource, 
            PhaseAgentManager agentManager, AtomsetIteratorMolecule neighborIterator) {
        super(space, phaseAgentSource, agentManager);
        singletAtomIterator = new AtomIteratorSinglet();
		singletPairIterator = new AtomsetIteratorSinglet(2);
        this.neighborIterator = neighborIterator;
	}
    
    public void setSimulation(Simulation sim) {
        phaseAgentManager.setRoot(sim.speciesRoot);
    }
    
    /**
     * @return Returns the cellRange.
     */
    public int getCellRange() {
        return cellRange;
    }

    /**
     * @param cellRange The cellRange to set.
     */
    public void setCellRange(int cellRange) {
        this.cellRange = cellRange;
    }

    protected void addToPotentialTypeList(PotentialAtomTypeWrapper wrapper) {
        boolean found = false;
        NeighborCriterion criterion = wrapper.potential.getCriterion();
        for (int i=0; i<criteriaArray.length; i++) {
            if (criteriaArray[i] == criterion) {
                found = true;
                break;
            }
        }
        if (!found) {
            criteriaArray = (NeighborCriterion[]) Arrays.addObject(criteriaArray, criterion);
        }
        super.addToPotentialTypeList(wrapper);
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
        if (!enabled)
            return;
        for (int i=0; i<criteriaArray.length; i++) {
            criteriaArray[i].setPhase(phase);
        }
        currentCellManager = (PhaseCellManager)phaseAgentManager.getAgents()[phase.getIndex()];
        AtomSet targetAtoms = id.getTargetAtoms();
        if (targetAtoms.count() == 0) {
            //no target atoms specified -- do one-target algorithm to
            // SpeciesMaster
            neighborIterator.setPhase(phase);
            neighborIterator.setDirection(IteratorDirective.Direction.UP);
            calculate(phase.getSpeciesMaster(), idUp, pc, getPotentials(
                    phase.getSpeciesMaster().type).getPotentials());
            if (lrcMaster != null) {
                lrcMaster.calculate(phase, id, pc);
            }
        } else if (targetAtoms instanceof Atom) {
            // one target atom
            neighborIterator.setPhase(phase);
            neighborIterator.setDirection(id.direction());
            calculate((Atom) targetAtoms, id, pc, getPotentials(
                    ((Atom) targetAtoms).type).getPotentials());
            if (lrcMaster != null) {
                lrcMaster.calculate(phase, id, pc);
            }
        } else {
            //more than one target atom
            super.calculate(phase, id, pc);
        }
    }//end calculate
	
    /**
     * Performs given PotentialCalculation using potentials/neighbors associated
     * with the given atom (if any).  Then, if atom is not a leaf atom, iteration over
     * child atoms is performed and process is repeated (recursively) with each on down
     * the hierarchy until leaf atoms are reached.
     */
    //TODO make a "TerminalGroup" indicator in type that permits child atoms but indicates that no potentials apply directly to them
	protected void calculate(Atom atom, IteratorDirective id, PotentialCalculation pc, final Potential[] potentials) {
        prepNbrIterator(atom);
        for(int i=0; i<potentials.length; i++) {
            switch (potentials[i].nBody()) {
            case 1:
                singletAtomIterator.setAtom(atom);
                pc.doCalculation(singletAtomIterator, id, potentials[i]);
                break;
            case 2:
                Potential2 p2 = (Potential2) potentials[i];
                NeighborCriterion nbrCriterion = p2.getCriterion();
                neighborIterator.setTarget(atom);
                neighborIterator.reset();
                while (neighborIterator.hasNext()) {
                    AtomPair pair = (AtomPair)neighborIterator.next();
                    if (nbrCriterion.accept(pair)) {
                        singletPairIterator.setAtom(pair);
                        pc.doCalculation(singletPairIterator, id, p2);
                    }
                }
                break;
            }
        }
            
		//if atom has children, repeat process with them
		if(!atom.node.isLeaf()) {
            //cannot use AtomIterator field because of recursive call
            AtomArrayList list = ((AtomTreeNodeGroup) atom.node).childList;
            int size = list.size();
            for (int i=0; i<size; i++) {
                Atom a = list.get(i);
                Potential[] childPotentials = getPotentials(a.type).getPotentials();
                calculate(a, id, pc, childPotentials);//recursive call
            }
		}
	}
    
    protected void prepNbrIterator(Atom atom) {
        ((Api1ASite)neighborIterator).setCentralSite(((NeighborSiteManager)currentCellManager).getSite(atom));
    }
    
    private final AtomIteratorSinglet singletAtomIterator;
	private final AtomsetIteratorSinglet singletPairIterator;
    private int cellRange;
    private final IteratorDirective idUp = new IteratorDirective();
    protected final AtomsetIteratorMolecule neighborIterator;
    protected PhaseCellManager currentCellManager;
    private NeighborCriterion[] criteriaArray = new NeighborCriterion[0];
    
    public static class PhaseAgentSiteManager implements PhaseAgentSource {
        public PhaseAgentSiteManager(int nCells) {
            this.nCells = nCells;
        }
        
        public Class getAgentClass() {
            return NeighborSiteManager.class;
        }
        
        public Object makeAgent(Phase phase) {
            return new NeighborSiteManager(phase,nCells);
        }
        
        public void releaseAgent(Object agent) {
        }
        
        private final int nCells;
    }
}
