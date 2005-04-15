/*
 * History
 * Created on Sep 20, 2004 by kofke
 */
package etomica.nbratom.cell;

import etomica.Atom;
import etomica.AtomPair;
import etomica.AtomSet;
import etomica.IteratorDirective;
import etomica.Phase;
import etomica.Potential;
import etomica.PotentialMaster;
import etomica.Simulation;
import etomica.Space;
import etomica.atom.AtomLinker;
import etomica.atom.AtomList;
import etomica.atom.AtomPositionDefinition;
import etomica.atom.AtomPositionDefinitionSimple;
import etomica.atom.AtomSequencerFactory;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.iterator.AtomsetIteratorSinglet;
import etomica.potential.Potential2;
import etomica.potential.PotentialCalculation;
import etomica.utility.Arrays;

/**
 * PotentialMaster used to implement neighbor listing.  Instance of this
 * class is given as an argument to the Simulation constructor.
 * Criteria specifying whether two atoms are neighbors for a particular potential
 * are specified in the setSpecies method of this class.
 * <br>
 */
public class PotentialMasterCell extends PotentialMaster {

	/**
	 * Invokes superclass constructor, specifying IteratorFactoryCell
     * for generating molecule iterators.  Sets default nCells of 10. 
	 */
	public PotentialMasterCell(Space space) {
        super(space,new IteratorFactoryCell());
        setNCells(10);
		singletIterator = new AtomsetIteratorSinglet(2);
        positionDefinition = new AtomPositionDefinitionSimple();
        neighborIterator = new Api1ACell(space.D());
	}
    
    /**
     * Performs cell-assignment potentialCalculation.  Assigns all molecules
     * to their cells, and invokes superclass method causing setup to be
     * performed iterating using species/potential hierarchy.
     */
    public void calculate(Phase phase, PotentialCalculationAgents pc) {
        super.calculate(phase, new IteratorDirective(), pc);
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
    	AtomSet targetAtoms = id.getTargetAtoms();
    	if (targetAtoms == null) {
    		//no target atoms specified -- do one-target algorithm to SpeciesMaster
            neighborIterator.setPhase(phase);
            neighborIterator.setDirection(IteratorDirective.UP);
    		calculate(phase.speciesMaster, idUp, pc, phase.speciesMaster.type.getNbrManagerAgent().getPotentials());
            if(lrcMaster != null) {
                lrcMaster.calculate(phase, id, pc);
            }
   	}
    	else if (targetAtoms instanceof Atom) {
    		// one target atom
            neighborIterator.setPhase(phase);
            neighborIterator.setDirection(id.direction());
			calculate((Atom)targetAtoms, id, pc, ((Atom)targetAtoms).type.getNbrManagerAgent().getPotentials());
            if(lrcMaster != null) {
                lrcMaster.calculate(phase, id, pc);
            }
    	}
    	else {
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
    //TODO make a "TerminalGroup" node that permits child atoms but indicates that no potentials apply directly to them
	private void calculate(Atom atom, IteratorDirective id, PotentialCalculation pc, Potential[] potentials) {
        int length = potentials.length;
		if (length > 0) {
            neighborIterator.setTarget(atom);
            neighborIterator.reset();
            while (neighborIterator.hasNext()) {
                AtomPair pair = neighborIterator.nextPair();
                singletIterator.setAtom(pair);
                for (int i=0; i<potentials.length; i++) {
                    if (((Potential2)potentials[i]).getCriterion().accept(pair)) {
                        pc.doCalculation(singletIterator,id,potentials[i]);
                    }
                }
            }
		}
		//if atom has children, repeat process with them
		if(!atom.node.isLeaf()) {
            //cannot use AtomIterator field because of recursive call
            AtomList list = ((AtomTreeNodeGroup) atom.node).childList;
            AtomLinker link = list.header.next;
            if (link != list.header) {
                Potential[] childPotentials = link.atom.type.getNbrManagerAgent().getPotentials();
                for (link = list.header.next; link != list.header; link = link.next) {
                    calculate(link.atom, id, pc, childPotentials);//recursive call
                }
            }
		}
	}
    
    public NeighborCellManager getNbrCellManager(Phase phase) {
        if(phase.getIndex() > neighborCellManager.length-1) {
            neighborCellManager = (NeighborCellManager[])Arrays.resizeArray(neighborCellManager, phase.getIndex()+1);
        }
        if(neighborCellManager[phase.getIndex()] == null) {
            neighborCellManager[phase.getIndex()] = new NeighborCellManager(phase,nCells,positionDefinition);
        }
        return neighborCellManager[phase.getIndex()];
    }

    public int getNCells() {
        return nCells;
    }
    public void setNCells(int cells) {
        nCells = cells;
    }
    
    public AtomSequencerFactory sequencerFactory() {return AtomSequencerCell.FACTORY;}
    
    public void setAtomPositionDefinition(AtomPositionDefinition positionDefinition) {
        this.positionDefinition = positionDefinition;
    }
    public AtomPositionDefinition getAtomPositionDefinition() {
        return positionDefinition;
    }

    public void setRange(double d) {
        neighborIterator.getNbrCellIterator().setRange(d);
    }
    
	private final AtomsetIteratorSinglet singletIterator;
    private NeighborCellManager[] neighborCellManager = new NeighborCellManager[0];
    private int nCells;
    private final IteratorDirective idUp = new IteratorDirective();
    private AtomPositionDefinition positionDefinition;
    private final Api1ACell neighborIterator;
    
}
