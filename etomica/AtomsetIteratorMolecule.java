/*
 * History
 * Created on Aug 31, 2004 by kofke
 */
package etomica;

/**
 * @author kofke
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public final class AtomsetIteratorMolecule extends AtomsetIteratorAdapter
		implements AtomsetIteratorPhaseDependent, AtomsetIteratorTargetDependent {

	private final Species[] species;
	private final AtomsetIteratorBasisDependent basisIterator;
	private final SpeciesAgent[] agents;
	private final int basisSize;
	
	/**
	 * @param iterator
	 */
	public AtomsetIteratorMolecule(Species[] species) {
		super(makeIterator(species));
		this.species = (Species[])species.clone();
		basisIterator = (AtomsetIteratorBasisDependent)iterator;
		basisSize = species.length;
		agents = new SpeciesAgent[basisSize];
	}

    public void setPhase(Phase phase) {
    	if (phase == null) {
    		basisIterator.setBasis(null);
    	}
    	else {
    		for (int i=0; i<basisSize; i++) {
    			agents[i] = phase.getAgent(species[i]);
    		}
    		basisIterator.setBasis(agents);
    	}
    }
    
    public void setTarget(Atom[] targetAtoms) {
    	basisIterator.setTarget(targetAtoms);
    }
    
    private static AtomsetIterator makeIterator(Species[] species) {
    	if (species == null || species.length == 0) {
    		throw new IllegalArgumentException("invalid species");
    	}
    	if (species.length==1) {
    		return new AtomIteratorDirectable();
    	}
    	if (species[0] == species[1]) {
    		return new ApiIntragroup();
    	}
    	return new ApiIntergroup();
    }
}
