/*
 * History
 * Created on Sep 20, 2004 by kofke
 */
package etomica.nbr;

import java.util.ArrayList;

import etomica.ApiInnerFixed;
import etomica.ApiMolecule;
import etomica.Atom;
import etomica.AtomArrayList;
import etomica.AtomIteratorArrayList;
import etomica.AtomIteratorListSimple;
import etomica.AtomIteratorSinglet;
import etomica.AtomSequencerFactory;
import etomica.AtomTreeNodeGroup;
import etomica.AtomsetIteratorMolecule;
import etomica.Default;
import etomica.IteratorDirective;
import etomica.NearestImageTransformerVector;
import etomica.Phase;
import etomica.Potential;
import etomica.Potential2;
import etomica.PotentialCalculation;
import etomica.PotentialMaster;
import etomica.Simulation;
import etomica.Space;
import etomica.Species;
import etomica.nbr.cell.AtomsetIteratorCellular;
import etomica.nbr.cell.IteratorFactoryCell;
import etomica.nbr.cell.NeighborCellManager;
import etomica.utility.Arrays;

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
        setMaxNeighborRange(Default.POTENTIAL_CUTOFF_FACTOR*1.5);
		neighborManager = new NeighborManager(this);
		atomIterator = new MyIterator();
		singletIterator = new AtomIteratorSinglet();
		pairIterator = new ApiInnerFixed(singletIterator, atomIterator);
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
    	Atom[] targetAtoms = id.getTargetAtoms();
    	if (targetAtoms.length == 0 || targetAtoms[0] == null) {
    		//no target atoms specified -- do one-target algorithm to SpeciesMaster
    		calculate(phase.speciesMaster, idUp, pc, phase.speciesMaster.type.getNbrManagerAgent().getPotentials());
    	}
    	else if (targetAtoms.length == 1 || targetAtoms[1] == null) {
    		// one target atom
			calculate(targetAtoms[0], id, pc, targetAtoms[0].type.getNbrManagerAgent().getPotentials());
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
			if (direction == IteratorDirective.UP || direction == null) {
				list = seq.getUpList();
                vectors = seq.getUpListNearestImageVector();
                atomIterator.nearestImageTransformer.setPlus(false);
//              list.length may be less than potentials.length, if atom hasn't yet interacted with another using one of the potentials
				for (int i=0; i<list.length; i++) {
					atomIterator.setList(list[i]);
                    atomIterator.setVectors(vectors[i]);
                    ((Potential2)potentials[i]).setNearestImageTransformer(atomIterator.nearestImageTransformer);
					//System.out.println("Up :"+atomIterator.size());
					pc.doCalculation(pairIterator, id, potentials[i]);
				}
			}
			if (direction == IteratorDirective.DOWN || direction == null) {
				list = seq.getDownList();
                vectors = seq.getDownListNearestImageVector();
                atomIterator.nearestImageTransformer.setPlus(true);
				for (int i=0; i<list.length; i++) {
					atomIterator.setList(list[i]);
                    atomIterator.setVectors(vectors[i]);
                    ((Potential2)potentials[i]).setNearestImageTransformer(atomIterator.nearestImageTransformer);
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
            if (listIterator.hasNext()) {
                Potential[] childPotentials = listIterator.peek()[0].type.getNbrManagerAgent().getPotentials();
                while(listIterator.hasNext()) {
                    calculate(listIterator.nextAtom(), id, pc, childPotentials);//recursive call
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
    public void setSpecies(Potential potential, Species[] species) {
    	if(potential.nBody() == 2) {
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
        ApiMolecule iterator = (ApiMolecule)iteratorFactory.makeMoleculeIterator(species);
        ((AtomsetIteratorCellular)iterator.getApiAA()).getNbrCellIterator().setRange(maxNeighborRange);
        ((AtomsetIteratorCellular)iterator.getApi1A()).getNbrCellIterator().setRange(maxNeighborRange);
        AtomsetIteratorMolecule iteratorFiltered = new ApiFiltered(iterator, (NeighborCriterionSimple)criterion);
		neighborManager.addCriterion(criterion);//add criterion to manager so criterion can be informed of the phase
    	for(int i=0; i<species.length; i++) {
    		species[i].moleculeFactory().getType().getNbrManagerAgent().addCriterion(criterion);//addCriterion method will prevent multiple additions of same criterion, if species are same
    	}
    	addPotential(potential, iteratorFiltered);
    }
    
    public void setSimulation(Simulation sim) {
//        sim.elementCoordinator.addMediatorPair(new etomica.Mediator.IntegratorPhase.NoCentralImage(sim.elementCoordinator));
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
    
    public void setMaxNeighborRange(double r) {
        maxNeighborRange = r;
    }
    
    public AtomSequencerFactory sequencerFactory() {return AtomSequencerNbr.FACTORY;}

	private final MyIterator atomIterator;
	private final AtomIteratorSinglet singletIterator;
	private final ApiInnerFixed pairIterator;
	private final NeighborManager neighborManager;
    private NeighborCellManager[] neighborCellManager = new NeighborCellManager[0];
    private int nCells;
    private double maxNeighborRange;
    private final IteratorDirective idUp = new IteratorDirective();
    private ArrayList[] vectors;
    
    public static class MyIterator extends AtomIteratorArrayList {
        
        public Atom nextAtom() {
            nearestImageTransformer.setNearestImageVector((Space.Vector)vector.get(cursor));
            Atom atom = super.nextAtom();
            return atom;
        }
        
        //TODO allAtoms
        
        public void setVectors(ArrayList vector) {
            this.vector = vector;
        }
        
        ArrayList vector;
        NearestImageTransformerVector nearestImageTransformer = new NearestImageTransformerVector();
    }
}
