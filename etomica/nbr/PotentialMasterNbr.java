/*
 * History
 * Created on Sep 20, 2004 by kofke
 */
package etomica.nbr;

import etomica.Atom;
import etomica.AtomSet;
import etomica.IteratorDirective;
import etomica.Phase;
import etomica.Potential;
import etomica.PotentialGroup;
import etomica.PotentialMaster;
import etomica.Space;
import etomica.Species;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomLinker;
import etomica.atom.AtomList;
import etomica.atom.AtomPositionDefinition;
import etomica.atom.AtomPositionDefinitionSimple;
import etomica.atom.AtomSequencerFactory;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.iterator.ApiInnerFixed;
import etomica.atom.iterator.ApiMolecule;
import etomica.atom.iterator.AtomIteratorArrayList;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.atom.iterator.AtomsetIteratorMolecule;
import etomica.nbr.cell.AtomsetIteratorCellular;
import etomica.nbr.cell.IteratorFactoryCell;
import etomica.nbr.cell.NeighborCellManager;
import etomica.potential.PotentialCalculation;
import etomica.utility.Arrays;
import etomica.utility.ObjectArrayList;

/**
 * PotentialMaster used to implement neighbor listing.  Instance of this
 * class is given as an argument to the Simulation constructor.
 * Criteria specifying whether two atoms are neighbors for a particular potential
 * are specified in the setSpecies method of this class.
 * <br>
 */
public class PotentialMasterNbr extends PotentialMaster {

	/**
	 * Invokes superclass constructor, specifying IteratorFactoryCell
     * for generating molecule iterators.  Sets default nCells of 10. 
	 */
	public PotentialMasterNbr(Space space) {
        super(space,new IteratorFactoryCell(space.D()));
        setNCells(10);
		neighborManager = new NeighborManager(this);
		atomIterator = new AtomIteratorArrayList();
		singletIterator = new AtomIteratorSinglet();
		pairIterator = new ApiInnerFixed(singletIterator, atomIterator);
        positionDefinition = new AtomPositionDefinitionSimple();
	}

    /**
     * Performs neighbor-setup potentialCalculation.  Assigns all molecules
     * to their cells, and invokes superclass method causing setup to be
     * performed iterating using species/potential hierarchy.
     */
    public void calculate(Phase phase, IteratorDirective id, PotentialCalculationNbrSetup pc) {
        getNbrCellManager(phase).assignCellAll();
        super.calculate(phase, id, pc);
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
    		calculate(phase.speciesMaster, idUp, pc, phase.speciesMaster.type.getNbrManagerAgent().getPotentials());
            if(lrcMaster != null) {
                lrcMaster.calculate(phase, id, pc);
            }
   	}
    	else if (targetAtoms instanceof Atom) {
    		// one target atom
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
            AtomSequencerNbr seq = (AtomSequencerNbr)atom.seq;
			singletIterator.setAtom(atom);
			IteratorDirective.Direction direction = id.direction();
			AtomArrayList[] list;
            ObjectArrayList[] vectors = null;
			if (direction == IteratorDirective.UP || direction == null) {
				list = seq.getUpList();
//              list.length may be less than potentials.length, if atom hasn't yet interacted with another using one of the potentials
				for (int i=0; i<list.length; i++) {
					atomIterator.setList(list[i]);
					//System.out.println("Up :"+atomIterator.size());
					pc.doCalculation(pairIterator, id, potentials[i]);
				}
			}
			if (direction == IteratorDirective.DOWN || direction == null) {
				list = seq.getDownList();
				for (int i=0; i<list.length; i++) {
					atomIterator.setList(list[i]);
					//System.out.println("Dn :"+atomIterator.size());
					pc.doCalculation(pairIterator, id, potentials[i]);
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

	//TODO update nbr radius for all criteria as more are added
    /**
     * Identifies given potential to apply to the given set of species,
     * and (if potential is a 2-body potential) constructs a default 
     * NeighborCriterionSimple instance to define the neighbors.
     */
    //TODO warn that potentialGroups don't get cells or criterion through this method
    public void setSpecies(Potential potential, Species[] species) {
    	if(potential.nBody() == 2 && !(potential instanceof PotentialGroup)) {
		    NeighborCriterion criterion = new NeighborCriterionSimple(space,potential.getRange(),2.0*potential.getRange());
	    	setSpecies(potential, species, criterion);
    	} else {
    	   	super.setSpecies(potential, species);
    	}
    }
    
    /**
     * Identifies given potential to apply to the given set of species, and 
     * sets the NeighborCriterion instance to define the neighbors.
     */
    public void setSpecies(Potential potential, Species[] species, NeighborCriterion criterion) {
    	if (species.length <= 1 || potential.nBody() != species.length) {
    		throw new IllegalArgumentException("Illegal species length");
    	}
        double neighborRange = criterion.getNeighborRange();
        ApiMolecule iterator = (ApiMolecule)iteratorFactory.makeMoleculeIterator(species);
        ((AtomsetIteratorCellular)iterator.getApiAA()).getNbrCellIterator().setRange(neighborRange);
        ((AtomsetIteratorCellular)iterator.getApi1A()).getNbrCellIterator().setRange(neighborRange);
        AtomsetIteratorMolecule iteratorFiltered = new ApiFiltered(iterator, criterion);
		neighborManager.addCriterion(criterion);//add criterion to manager so criterion can be informed of the phase
    	for(int i=0; i<species.length; i++) {
            species[i].moleculeFactory().getType().getNbrManagerAgent().addCriterion(criterion);//addCriterion method will prevent multiple additions of same criterion, if species are same
        }
    	addPotential(potential, iteratorFiltered, moleculeTypes(species));
    }
    
    /**
     * Identifies given potential to apply to the given set of species, and 
     * sets the NeighborCriterion instance to define the neighbors.
     */
    public void setSpecies(Potential potential, Species[] species, double potentialRange) {
        if (species.length <= 1 || potential.nBody() != species.length) {
            throw new IllegalArgumentException("Illegal species length");
        }
        ApiMolecule iterator = (ApiMolecule)iteratorFactory.makeMoleculeIterator(species);
        ((AtomsetIteratorCellular)iterator.getApiAA()).getNbrCellIterator().setRange(potentialRange);
        ((AtomsetIteratorCellular)iterator.getApi1A()).getNbrCellIterator().setRange(potentialRange);
        addPotential(potential, iterator, moleculeTypes(species));
    }
    
    public NeighborManager getNeighborManager() {return neighborManager;}

    
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
    
    public AtomSequencerFactory sequencerFactory() {return AtomSequencerNbr.FACTORY;}
    
    public void setAtomPositionDefinition(AtomPositionDefinition positionDefinition) {
        this.positionDefinition = positionDefinition;
    }
    public AtomPositionDefinition getAtomPositionDefinition() {
        return positionDefinition;
    }

	private final AtomIteratorArrayList atomIterator;
	private final AtomIteratorSinglet singletIterator;
	private final ApiInnerFixed pairIterator;
	private final NeighborManager neighborManager;
    private NeighborCellManager[] neighborCellManager = new NeighborCellManager[0];
    private int nCells;
    private final IteratorDirective idUp = new IteratorDirective();
    private AtomPositionDefinition positionDefinition;
    
}
