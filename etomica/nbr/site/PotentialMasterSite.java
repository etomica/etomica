/*
 * History
 * Created on Sep 20, 2004 by kofke
 */
package etomica.nbr.site;

import etomica.IteratorDirective;
import etomica.Phase;
import etomica.Potential;
import etomica.PotentialMaster;
import etomica.Space;
import etomica.atom.Atom;
import etomica.atom.AtomLinker;
import etomica.atom.AtomList;
import etomica.atom.AtomPair;
import etomica.atom.AtomPositionDefinition;
import etomica.atom.AtomSequencerFactory;
import etomica.atom.AtomSet;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.atom.iterator.AtomsetIteratorMolecule;
import etomica.atom.iterator.AtomsetIteratorSinglet;
import etomica.nbr.NeighborCriterion;
import etomica.nbr.PotentialCalculationAgents;
import etomica.nbr.cell.IteratorFactoryCell;
import etomica.nbr.cell.NeighborCellManager;
import etomica.potential.Potential2;
import etomica.potential.PotentialCalculation;

/**
 * 
 * <br>
 */
public class PotentialMasterSite extends PotentialMaster {

	/**
	 * Invokes superclass constructor, specifying IteratorFactoryCell
     * for generating molecule iterators.  Sets default nCells of 10 and
     * position definition to null, so that atom type's definition is used
     * to assign cells. 
	 */
	public PotentialMasterSite(Space space) {
        this(space, null);
    }
    
    /**
     * Constructs class using given position definition for all atom cell assignments.
     * @param positionDefinition if null, specifies use of atom type's position definition
     */
    public PotentialMasterSite(Space space, AtomPositionDefinition positionDefinition) {
        this(space,positionDefinition,new Api1ASite(space.D()));
    }
    
    public PotentialMasterSite(Space space, AtomPositionDefinition positionDefinition, AtomsetIteratorMolecule neighborIterator) {
        super(space,new IteratorFactoryCell());
        setNCells(10);
        singletAtomIterator = new AtomIteratorSinglet();
		singletPairIterator = new AtomsetIteratorSinglet(2);
        this.positionDefinition = positionDefinition;//use atom type's position definition as default
        this.neighborIterator = neighborIterator;
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
        if (!enabled)
            return;
        AtomSet targetAtoms = id.getTargetAtoms();
        if (targetAtoms.count() == 0) {
            //no target atoms specified -- do one-target algorithm to
            // SpeciesMaster
            neighborIterator.setPhase(phase);
            neighborIterator.setDirection(IteratorDirective.UP);
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
	private void calculate(Atom atom, IteratorDirective id, PotentialCalculation pc, final Potential[] potentials) {
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
            
//        if (length > 0) {
//            neighborIterator.setTarget(atom);
//            neighborIterator.reset();
//            while (neighborIterator.hasNext()) {
//                AtomPair pair = neighborIterator.nextPair();
//                singletPairIterator.setAtom(pair);
//                for (int i=0; i<potentials.length; i++) {
//                    if (((Potential2)potentials[i]).getCriterion().accept(pair)) {
//                        pc.doCalculation(singletPairIterator,id,potentials[i]);
//                    }
//                }
//            }
//		}
		//if atom has children, repeat process with them
		if(!atom.node.isLeaf()) {
            //cannot use AtomIterator field because of recursive call
            AtomList list = ((AtomTreeNodeGroup) atom.node).childList;
            AtomLinker link = list.header.next;
            if (link != list.header) {
                for (link = list.header.next; link != list.header; link = link.next) {
                    Potential[] childPotentials = getPotentials(link.atom.type).getPotentials();
                    calculate(link.atom, id, pc, childPotentials);//recursive call
                }
            }
		}
	}
    
    public NeighborCellManager getNbrCellManager(Phase phase) {
        NeighborCellManager manager = (NeighborCellManager)phase.getCellManager();
        if (manager == null) {
            manager = new NeighborCellManager(phase,nCells,positionDefinition);
            phase.setCellManager(manager);
        }
        return manager;
    }

    public int getNCells() {
        return nCells;
    }
    public void setNCells(int cells) {
        if (cells < 3) {
            throw new IllegalArgumentException("You must have at least 3 cells");
        }
        nCells = cells;
    }
    
    public AtomSequencerFactory sequencerFactory() {return AtomSequencerSite.FACTORY;}
    
    private final AtomIteratorSinglet singletAtomIterator;
	private final AtomsetIteratorSinglet singletPairIterator;
    private int nCells;
    private final IteratorDirective idUp = new IteratorDirective();
    private final AtomPositionDefinition positionDefinition;
    protected final AtomsetIteratorMolecule neighborIterator;
    
}
