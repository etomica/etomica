/*
 * History
 * Created on Sep 20, 2004 by kofke
 */
package etomica.nbr.list;

import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomLinker;
import etomica.atom.AtomList;
import etomica.atom.AtomPositionDefinition;
import etomica.atom.AtomSequencerFactory;
import etomica.atom.AtomSet;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.iterator.ApiInnerFixed;
import etomica.atom.iterator.AtomIteratorArrayList;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.atom.iterator.IteratorDirective;
import etomica.nbr.PotentialCalculationUpdateTypeList;
import etomica.nbr.cell.IteratorFactoryCell;
import etomica.nbr.cell.NeighborCellManager;
import etomica.phase.Phase;
import etomica.potential.Potential;
import etomica.potential.PotentialCalculation;
import etomica.potential.PotentialMaster;
import etomica.space.Space;

/**
 * PotentialMaster used to implement neighbor listing.  Instance of this
 * class is given as an argument to the Simulation constructor.
 * Criteria specifying whether two atoms are neighbors for a particular potential
 * are specified in the setSpecies method of this class.
 * <br>
 */
public class PotentialMasterNbr extends PotentialMaster {

	/**
     * Default constructor.  Range should be set manually via 
     * neighborCellManager.
     */
    public PotentialMasterNbr(Space space) {
        this(space,1.0);
    }
    
    /**
	 * Invokes superclass constructor, specifying IteratorFactoryCell
     * for generating molecule iterators.  Sets default nCells of 10. 
	 */
	public PotentialMasterNbr(Space space, double range) {
        this(space, range, null);
    }
    
    /**
     * Constructs class using given position definition for all atom cell assignments.
     * @param positionDefinition if null, specifies use of atom type's position definition
     */
    public PotentialMasterNbr(Space space, double range, AtomPositionDefinition positionDefinition) {
        super(space,new IteratorFactoryCell());
        neighborManager = new NeighborListManager(this, range);
        atomIterator = new AtomIteratorArrayList();
        singletIterator = new AtomIteratorSinglet();
        pairIterator = new ApiInnerFixed(singletIterator, atomIterator);
        swappedPairIterator = new ApiInnerFixed(singletIterator, atomIterator, true);
        this.positionDefinition = positionDefinition;
        cellRange = 2;
	}
    
    public void setRange(double range) {
        neighborManager.setRange(range);
    }
    
    public void getRange() {
        neighborManager.getRange();
    }
    
    /**
     * Performs cell-assignment potentialCalculation.  Assigns all molecules
     * to their cells, and invokes superclass method causing setup to be
     * performed iterating using species/potential hierarchy.
     */
    public void calculate(Phase phase, PotentialCalculationUpdateTypeList pc) {
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
    	if (targetAtoms.count() == 0) {
    		//no target atoms specified -- do one-target algorithm to SpeciesMaster
    		calculate(phase.getSpeciesMaster(), idUp, pc, getPotentials(phase.getSpeciesMaster().type).getPotentials());
            if(lrcMaster != null) {
                lrcMaster.calculate(phase, id, pc);
            }
    	}
    	else if (targetAtoms instanceof Atom) {
    		// one target atom
			calculate((Atom)targetAtoms, id, pc, getPotentials(((Atom)targetAtoms).type).getPotentials());
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
        singletIterator.setAtom(atom);
        IteratorDirective.Direction direction = id.direction();
        for(int i=0; i<potentials.length; i++) {
            switch (potentials[i].nBody()) {
            case 1:
                pc.doCalculation(singletIterator, id, potentials[i]);
                break;
            case 2:
                AtomSequencerNbr seq = (AtomSequencerNbr)atom.seq;
                AtomArrayList[] list;
                if (direction != IteratorDirective.DOWN) {
                    list = seq.getUpList();
//                  list.length may be less than potentials.length, if atom hasn't yet interacted with another using one of the potentials
                    if(i < list.length) {
                        atomIterator.setList(list[i]);
                        //System.out.println("Up :"+atomIterator.size());
                        pc.doCalculation(pairIterator, id, potentials[i]);
                    }
                }
                if (direction != IteratorDirective.UP) {
                    list = seq.getDownList();
                    if(i < list.length) {
                        atomIterator.setList(list[i]);
                        //System.out.println("Dn :"+atomIterator.size());
                        pc.doCalculation(swappedPairIterator, id, potentials[i]);
                    }
                }
                break;//switch
            }//end of switch
        }//end of for
        
//        if (length > 0) {
//            AtomSequencerNbr seq = (AtomSequencerNbr)atom.seq;
//			singletIterator.setAtom(atom);
//			IteratorDirective.Direction direction = id.direction();
//			AtomArrayList[] list;
//			if (direction == IteratorDirective.UP || direction == null) {
//				list = seq.getUpList();
////              list.length may be less than potentials.length, if atom hasn't yet interacted with another using one of the potentials
//				for (int i=0; i<list.length; i++) {
//					atomIterator.setList(list[i]);
//					//System.out.println("Up :"+atomIterator.size());
//					pc.doCalculation(pairIterator, id, potentials[i]);
//				}
//			}
//			if (direction == IteratorDirective.DOWN || direction == null) {
//				list = seq.getDownList();
//				for (int i=0; i<list.length; i++) {
//					atomIterator.setList(list[i]);
//					//System.out.println("Dn :"+atomIterator.size());
//					pc.doCalculation(pairIterator, id, potentials[i]);
//				}
//			}
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

    public NeighborListManager getNeighborManager() {return neighborManager;}

    
    public NeighborCellManager getNbrCellManager(Phase phase) {
        NeighborCellManager manager = (NeighborCellManager)phase.getCellManager();
        if (manager == null) {
            manager = new NeighborCellManager(phase, neighborManager.getRange(), positionDefinition);
            phase.setCellManager(manager);
        }
        else {
            manager.setPotentialRange(neighborManager.getRange());
        }
        manager.setCellRange(cellRange);
        return manager;
    }

    public void setCellRange(int newCellRange) {
        cellRange = newCellRange;
    }

    public int getCellRange() {
        return cellRange;
    }
    
    public AtomSequencerFactory sequencerFactory() {return AtomSequencerNbr.FACTORY;}
    
    public AtomPositionDefinition getAtomPositionDefinition() {
        return positionDefinition;
    }

	private final AtomIteratorArrayList atomIterator;
	private final AtomIteratorSinglet singletIterator;
	private final ApiInnerFixed pairIterator;
    private final ApiInnerFixed swappedPairIterator;
	private final NeighborListManager neighborManager;
    private int cellRange;
    private final IteratorDirective idUp = new IteratorDirective();
    private final AtomPositionDefinition positionDefinition;
    
}
