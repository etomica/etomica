/*
 * History
 * Created on Aug 31, 2004 by kofke and schultz
 */
package etomica.atom.iterator;

import etomica.atom.Atom;
import etomica.atom.AtomsetArray;
import etomica.atom.SpeciesAgent;
import etomica.phase.Phase;
import etomica.species.Species;

/**
 * Singlet iterator that returns an atom set formed from one or more of the
 * species agents of a phase. Species must be specified at construction, and
 * when a subsequent call to setPhase is made, the basis for iteration will be
 * set to the species agents in that phase.  This iterator is used for calculations
 * involving zero-body potentials.
 */
public final class AtomsetIteratorSpeciesAgent extends AtomsetIteratorAdapter
		implements AtomsetIteratorPDT {

	/**
	 * Species are specified at construction and cannot be changed afterward.
     * Size of atomset is equal to the number of species given here.
	 */
	public AtomsetIteratorSpeciesAgent(Species[] species) {
		super(new AtomsetIteratorSinglet(species.length));
		this.species = (Species[])species.clone();
		basisSize = species.length;
		agents = new SpeciesAgent[basisSize];
        agentUsed = new boolean[basisSize];
        atoms = new AtomsetArray(agents);
        ((AtomsetIteratorSinglet)iterator).setAtom(atoms);
	}

	/**
     * Specifies the phase for which the pre-specified species agents are given.
     * If given phase is null, no iterates will be returned until a phase is
     * given via another call to this method.
     */
    public void setPhase(Phase phase) {
        if(this.phase != phase) {
            this.phase = phase;
            if (phase != null) {
                for (int i=0; i<basisSize; i++) {
                    agents[i] = phase.getAgent(species[i]);
                }
            }
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

    /**
     * Allows iteration only if each atom in given AtomSet is
     * derived from one of the species, with each species being
     * applied to no more than one target atom.
     */
    public void setTarget(Atom newTargetAtom) {
        if (newTargetAtom == null) {
            canReset = true;
        }
        else {
            canReset = true;
        	for(int j=0; j<basisSize; j++) {
                agentUsed[j] = false;
            }
    		boolean match = false;
    		for (int j=0; j<basisSize; j++) {
    			if (!agentUsed[j] && newTargetAtom.getNode().isDescendedFrom(agents[j])) {
    				agentUsed[j] = true;
    				match = true;
    				break;
    			}
    		}
    		if (!match) {
    			canReset = false;
    		}
        }
    }
    
    /**
     * Has no effect.
     */
    public void setDirection(IteratorDirective.Direction direction) { }

    private static final long serialVersionUID = 1L;
    private final Species[] species;
	private final SpeciesAgent[] agents;
	private final int basisSize;
	private boolean canReset = true;
    private final boolean[] agentUsed;
	private Phase phase;
    private final AtomsetArray atoms;
	
}
