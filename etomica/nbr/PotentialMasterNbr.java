/*
 * History
 * Created on Sep 20, 2004 by kofke
 */
package etomica.nbr;

import etomica.ApiInnerFixed;
import etomica.Atom;
import etomica.AtomArrayList;
import etomica.AtomIteratorArrayList;
import etomica.AtomIteratorListSimple;
import etomica.AtomIteratorSinglet;
import etomica.AtomSequencerFactory;
import etomica.AtomTreeNodeGroup;
import etomica.AtomsetIterator;
import etomica.AtomsetIteratorFiltered;
import etomica.AtomsetIteratorMolecule;
import etomica.IteratorDirective;
import etomica.Phase;
import etomica.Potential;
import etomica.PotentialCalculation;
import etomica.PotentialGroup;
import etomica.PotentialMaster;
import etomica.Simulation;
import etomica.Space;
import etomica.Species;
import etomica.nbr.cell.IteratorFactoryCell;
import etomica.nbr.cell.NeighborCellManager;
import etomica.utility.Arrays;

/**
 * PotentialMaster used to implement neighbor listing.  Instance of this
 * class is given as an argument to the Simulation constructor.
 * Criteria specifying whether two atoms are neighbors for a particular potential
 * are specified in the setSpecies method of this class.
 * <br>
 * Neighbor-list facility is implemented as follows.  Neighbor lists are held
 * by the atom's sequencer (AtomSequencerNbr).  Lists are keyed to the potential,
 * so given an atom it is possible to iterate over all current neighbors interacting
 * with it via a particular potential.  Such lists are kept only for "real" potentials,
 * not potential groups.  PotentialMaster constructs a NeighborManager, which listens
 * for interval events, and must therefore be registered with all integrators as
 * an interval listener. The neighborManager is responsible for keeping the neighbor
 * lists up to date.  The calculate method of PotentialMasterNbr is configured to 
 * perform neighbor-list iteration, or to update the neighbor lists if invoked with
 * a PotentialCalculationNbrSetup instance.
 */
public class PotentialMasterNbr extends PotentialMaster {

	/**
	 * @param space 
	 */
	public PotentialMasterNbr(Space space) {
        super(space,IteratorFactoryCell.INSTANCE);
        setNCells(10);
		neighborManager = new NeighborManager(this);
		atomIterator = new AtomIteratorArrayList();
		singletIterator = new AtomIteratorSinglet();
		pairIterator = new ApiInnerFixed(singletIterator, atomIterator);
	}

    /**
     * Overrides superclass method to enable direct neighbor-list iteration instead
     * of iteration via species/potential hierarchy.  If given PotentialCalculation
     * is for neighbor setup (instance of PotentialCalculationNbrSetup) then superclass
     * method is invoked and set-up is performed iterating using species/potential 
     * hierarchy.  Otherwise neighbor list facility is used to perform given calculation.
     * If no target atoms are specified in directive, neighborlist iteration is begun
     * with speciesMaster of phase, and repeated recursively down species hierarchy; if
     * one atom is specified, neighborlist iteration is performed on it and down
     * species hierarchy from it; if two or more atoms are specified, superclass method
     * is invoked.
     */
	public void calculate(Phase phase, IteratorDirective id, PotentialCalculation pc) {
		if(!enabled) return;
	   	if(pc instanceof PotentialCalculationNbrSetup) {
            getNbrCellManager(phase).assignCellAll();
	   		super.calculate(phase, id, pc);
	   	}
 		else {
 	    	Atom[] targetAtoms = id.getTargetAtoms();
 	    	if (targetAtoms.length == 0 || targetAtoms[0] == null) {
 	    		//no target atoms specified -- do one-target algorithm to SpeciesMaster
 	    		calculate(phase.speciesMaster, new IteratorDirective(), pc);
 	    	}
 	    	else if (targetAtoms.length == 1 || targetAtoms[1] == null) {
 	    		// one target atom
    			calculate(targetAtoms[0], id, pc);
 	    	}
 	    	else {
 	    		//more than one target atom
 	    		super.calculate(phase, id, pc);
 	    	}
		}
		
	}//end calculate
	
    /**
     * Performs given PotentialCalculation using potentials/neighbors associated
     * with the given atom (if any).  Then, if atom is not a leaf atom, iteration over
     * child atoms is performed and process is repeated (recursively) with each on down
     * the hierarchy until leaf atoms are reached.
     */
    //TODO make a "TerminalGroup" node that permits child atoms but indicates that no potentials apply directly to them
	private void calculate(Atom atom, IteratorDirective id, PotentialCalculation pc) {
		Potential[] potentials = atom.type.getNbrManagerAgent().getPotentials();
		int length = potentials.length;
		if (length > 0) {
            AtomSequencerNbr seq = (AtomSequencerNbr)atom.seq;
			singletIterator.setAtom(atom);
			IteratorDirective.Direction direction = id.direction();
			AtomArrayList[] list;
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
			//TODO if instantiation is expensive, try doing iteration explicitly (cannot use class variable because of recursive call)
			AtomIteratorListSimple listIterator = new AtomIteratorListSimple();
			listIterator.setList(((AtomTreeNodeGroup)atom.node).childList);
			listIterator.reset();
			while(listIterator.hasNext()) {
				calculate(listIterator.nextAtom(), id, pc);//recursive call
			}
		}
	}

	public void addPotentialNotify(Potential potential) {
		if(potential instanceof PotentialGroup) return;
		else {
			//TODO add to list of concrete potentials
		}
	}

    public void removePotentialNotify(Potential potential) {
		if(potential instanceof PotentialGroup) return;
		else {
			//TODO remove from list of concrete potentials
		}
	}
	

	//TODO update nbr radius for all criteria as more are added
    /**
     * Identifies given potential to apply to the given set of species,
     * and (if potential is a 2-body potential) constructs a default 
     * NeighborCriterionSimple instance to define the neighbors.
     */
    public void setSpecies(Potential potential, Species[] species) {
    	if(potential.nBody() == 2) {
		    NeighborCriterion criterion = new NeighborCriterionSimple(space,potential.getRange(),2.0*potential.getRange());
	    	setSpecies(potential, species, criterion);
    	} else {
    	   	AtomsetIterator iterator = new AtomsetIteratorMolecule(species,iteratorFactory);
    	    addPotential(potential, iterator);
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
    	AtomsetIterator iterator = new AtomsetIteratorMolecule(species,iteratorFactory);
    	iterator = new AtomsetIteratorFiltered(iterator, criterion);
		neighborManager.addCriterion(criterion);
    	for(int i=0; i<species.length; i++) {
    		species[i].moleculeFactory().getType().getNbrManagerAgent().addCriterion(criterion);//addCriterion method will prevent multiple additions of same criterion, if species are same
    	}
    	addPotential(potential, iterator);
    }
    
    public void setSimulation(Simulation sim) {
        sim.elementCoordinator.addMediatorPair(new etomica.Mediator.IntegratorPhase.NoCentralImage(sim.elementCoordinator));
    }

    public NeighborManager getNeighborManager() {return neighborManager;}

    
    public NeighborCellManager getNbrCellManager(Phase phase) {
        if(phase.index > neighborCellManager.length-1) {
            neighborCellManager = (NeighborCellManager[])Arrays.resizeArray(neighborCellManager, phase.index+1);
        }
        if(neighborCellManager[phase.index] == null) {
            neighborCellManager[phase.index] = new NeighborCellManager(phase,nCells);
        }
        return neighborCellManager[phase.index];
    }

    public int getNCells() {
        return nCells;
    }
    public void setNCells(int cells) {
        nCells = cells;
    }
    public AtomSequencerFactory sequencerFactory() {return AtomSequencerNbr.FACTORY;}

	private final AtomIteratorArrayList atomIterator;
	private final AtomIteratorSinglet singletIterator;
	private final ApiInnerFixed pairIterator;
	private final NeighborManager neighborManager;
    private NeighborCellManager[] neighborCellManager = new NeighborCellManager[0];
    private int nCells;
}
