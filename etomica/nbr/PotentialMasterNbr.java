/*
 * History
 * Created on Sep 20, 2004 by kofke
 */
package etomica.nbr;

import etomica.*;

/**
 * @author kofke
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class PotentialMasterNbr extends PotentialMaster {

	/**
	 * @param space
	 */
	public PotentialMasterNbr(Space space) {
		super(space);
		neighborManager = new NeighborManager(this);
		atomIterator = new AtomIteratorArrayList();
		singletIterator = new AtomIteratorSinglet();
		pairIterator = new ApiInnerFixed(singletIterator, atomIterator);
	}

	public void calculate(Phase phase, IteratorDirective id, PotentialCalculation pc) {
		if(!enabled) return;
	   	if(pc instanceof PotentialCalculationNbrSetup) {
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
				for (int i=0; i<length; i++) {
					atomIterator.setList(list[i]);
					//System.out.println("Up :"+atomIterator.size());
					pc.doCalculation(pairIterator, id, potentials[i]);
				}
			}
			if (direction == IteratorDirective.DOWN || direction == null) {
				list = seq.getDownList();
				for (int i=0; i<length; i++) {
					atomIterator.setList(list[i]);
					//System.out.println("Dn :"+atomIterator.size());
					pc.doCalculation(pairIterator, id, potentials[i]);
				}
			}
		}
		//if atom has children, repeat process with them
		if(!atom.node.isLeaf()) {
			//TODO if instantiation is expensive, try doing iteration explicitly
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
    public void setSpecies(Potential potential, Species[] species) {
    	if(potential.nBody() == 2) {
		    NeighborCriterion criterion = new NeighborCriterionSimple(space,potential.getRange(),2.0*potential.getRange());
	    	setSpecies(potential, species, criterion);
    	} else {
    	   	AtomsetIterator iterator = new AtomsetIteratorMolecule(species);
    	    addPotential(potential, iterator);
    	}
    }
    
    public void setSpecies(Potential potential, Species[] species, NeighborCriterion criterion) {
    	if (species.length <= 1 || potential.nBody() != species.length) {
    		throw new IllegalArgumentException("Illegal species length");
    	}
    	AtomsetIterator iterator = new AtomsetIteratorMolecule(species);
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

    public AtomSequencer.Factory sequencerFactory() {return AtomSequencerNbr.FACTORY;}

	private final AtomIteratorArrayList atomIterator;
	private final AtomIteratorSinglet singletIterator;
	private final ApiInnerFixed pairIterator;
	private final NeighborManager neighborManager;
}
