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
	 * @param nBody
	 */
	public PotentialMasterNbr(Simulation sim) {
		super(sim);
	}

	public void calculate(Phase phase, IteratorDirective id, PotentialCalculation pc) {
		if(!enabled) return;
	   	if(pc instanceof PotentialCalculationNbrSetup) {
	   		//TODO clear neighbor lists
	   		super.calculate(phase, id, pc);
	   	}
 		else {
 	    	Atom[] targetAtoms = id.getTargetAtoms();
 	    	switch(targetAtoms.length) {
 	    		case 0:
 	    	//no target atoms specified -- do one-target algorithm to SpeciesMaster
 	    			calculate(phase.speciesMaster, new IteratorDirective(), pc);
 	    			break;
 	    		case 1:
 	    	//one target atom specified -- 
 	    		//if it has potential/nbrs, loop through them
 	    		//if it has children, loop through them and do same
 	    			calculate(targetAtoms[0], id, pc);
 	    			break;
 	    		default:
 	    	//more than one target atom --
 	    		    super.calculate(phase, id, pc);
 	    			break;
 	    	}
		}
		
	}//end calculate
	
	private void calculate(Atom atom, IteratorDirective id, PotentialCalculation pc) {
		AtomSequencerNbr seq = (AtomSequencerNbr)atom.seq;
		Potential[] potentials = atom.type.getNbrManagerAgent().getPotentials();
		int length = potentials.length;
		if (length > 0) {
			IteratorDirective.Direction direction = id.direction();
			AtomArrayList[] list;
			if (direction == IteratorDirective.UP || direction == null) {
				list = seq.getUpList();
				for (int i=0; i<length; i++) {
					atomIterator.setList(list[i]);
					pc.doCalculation(atomIterator, id, potentials[i]);
				}
			}
			if (direction == IteratorDirective.DOWN || direction == null) {
				list = seq.getDownList();
				for (int i=0; i<length; i++) {
					atomIterator.setList(list[i]);
					pc.doCalculation(atomIterator, id, potentials[i]);
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
	
	private final AtomIteratorArrayList atomIterator = new AtomIteratorArrayList();
}
