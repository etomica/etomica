/*
 * History
 * Created on Aug 31, 2004 by kofke and schultz
 */
package etomica.atom.iterator;

import etomica.AtomSet;
import etomica.IteratorDirective;
import etomica.Phase;
import etomica.Species;
import etomica.SpeciesAgent;
import etomica.atom.AtomsetArray;

/**
 * Iterates molecules (children of SpeciesAgent atoms), as specified by phase
 * and species.  Species must be specified at construction, and when a 
 * subsequent call to setPhase is made, the basis for iteration will be 
 * set to the species agents in that phase.  
 */
//PotentialMaster uses instances of this class for all iterators it assigns
//to the potentials added to it via setSpecies.
public final class AtomsetIteratorSpeciesAgent extends AtomsetIteratorAdapter
		implements AtomsetIteratorMolecule {

	/**
	 * Constructor requires specification of species for which molecules
	 * are subject to iteration.  If a single species is given in the array,
	 * then iteration is performed over single-molecules of that species.  If
	 * two species are given in the array, then iteration is performed over pairs
	 * formed from molecules of the species; in this case the species may be the
	 * same instance (resulting in intra-species pair iteration) or different
	 * instances (resulting in inter-species pair iteration).  If null, zero, or more than
	 * two species are specified, an IllegalArgumentException is thrown.
	 */
	public AtomsetIteratorSpeciesAgent(Species[] species) {
		super(new AtomsetIteratorSinglet(species.length));
		this.species = (Species[])species.clone();
		basisSize = species.length;
		agents = new SpeciesAgent[basisSize];
        atoms = new AtomsetArray(basisSize);
	}

	/**
	 * Specifies the phase for which the pre-specified species' molecules
	 * will be subject to iteration.  If given phase is null, no iterates will
	 * be returned until a phase is given via another call to this method.
	 */
    public void setPhase(Phase phase) {
    	this.phase = phase;
    	if (phase != null) {
    		for (int i=0; i<basisSize; i++) {
    			agents[i] = phase.getAgent(species[i]);
    		}
            atoms.setAtoms(agents);
    		((AtomsetIteratorSinglet)iterator).setAtom(atoms);
    	}
    	unset();
    }

    public void reset() {
    	if (canReset && phase != null)
    		super.reset();
    	else
    		unset();
    }
    
    /**
     * Returns an array with the species that were set for this
     * iterator at its construction.  Copy of array is returned, so
     * any change to element of array will not affect iterator.
     */
    public Species[] getSpecies() {
    	return (Species[])species.clone();
    }

    public void setTarget(AtomSet targetAtoms) {
    	if (targetAtoms.count() > basisSize) {
    		canReset = false;
    		return;
    	}
    	boolean[] agentUsed = new boolean[basisSize];
    	for (int i=0; i<targetAtoms.count(); i++) {
    		boolean match = false;
    		for (int j=0; j<agents.length; j++) {
    			if (!agentUsed[j] && targetAtoms.getAtom(i).node.parentSpeciesAgent() == agents[j]) {
    				agentUsed[j] = true;
    				match = true;
    				break;
    			}
    		}
    		if (!match) {
    			canReset = false;
    			return;
    		}
    	}
    }
    
    public void setDirection(IteratorDirective.Direction direction) { }

    private final Species[] species;
	private final SpeciesAgent[] agents;
	private final int basisSize;
	private boolean canReset = true;
	private Phase phase;
    private final AtomsetArray atoms;
	
}
