/*
 * History
 * Created on Aug 31, 2004 by kofke and schultz
 */
package etomica;

/**
 * Iterates molecules (children of SpeciesAgent atoms), as specified by phase
 * and species.  Species must be specified at construction, and when a 
 * subsequent call to setPhase is made, the basis for iteration will be 
 * set to the species agents in that phase.  
 */
//PotentialMaster uses instances of this class for all iterators it assigns
//to the potentials added to it via setSpecies.
public final class AtomsetIteratorMolecule extends AtomsetIteratorAdapter
		implements AtomsetIteratorPhaseDependent, AtomsetIteratorTargetable,
		AtomsetIteratorDirectable {

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
	public AtomsetIteratorMolecule(Species[] species) {
		this(species, false);
	}
	public AtomsetIteratorMolecule(Species[] species, boolean dummy) {
		super(makeIterator(species));
		this.species = (Species[])species.clone();
		basisIterator = (AtomsetIteratorBasisDependent)iterator;
		basisSize = species.length;
		agents = new SpeciesAgent[basisSize];
		targetAtoms = new Atom[nBody()];
		
		//ignore direction specifications if this is a single-atom iterator
		//or if it is a intra-group pair iterator (in which case the direction
		//specification is passed on to the wrapped intra-group iterator).
		ignoreDirection = (basisSize == 1) || (species[0] == species[1]);
	}

	/**
	 * Specifies the phase for which the pre-specified species' molecules
	 * will be subject to iteration.  If given phase is null, no iterates will
	 * be returned until a phase is given via another call to this method.
	 */
    public void setPhase(Phase phase) {
    	if (phase == null) {
    		basisIterator.setBasis(null);
    		agents[0] = null;
    	}
    	else {
    		for (int i=0; i<basisSize; i++) {
    			agents[i] = phase.getAgent(species[i]);
    		}
    		needBasisUpdate = true;
    	}
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
    
    
    public void setTarget(Atom[] targetAtoms) {
    	System.arraycopy(targetAtoms,0,this.targetAtoms,0,targetAtoms.length);
    	for(int i=targetAtoms.length; i<this.targetAtoms.length; i++) {
    		this.targetAtoms[i] = null;
    	}
    	needBasisUpdate = true;
    	unset();
    }
    
    public void reset() {
    	if (needBasisUpdate && agents[0] != null) {
        	needBasisUpdate = false;
    		if (targetAtoms != null && targetAtoms[0] != null) {
    			SpeciesAgent targetAgent = targetAtoms[0].node.parentSpeciesAgent();
    			if (targetAgent != agents[0]) {
    				int i=basisSize;
    				if (!ignoreDirection) {
    					for (i=1; i<basisSize; i++) {
    						if (targetAgent == agents[i]) {
    							agents[i] = agents[0];
    							agents[0] = targetAgent;
    							break;
    						}
    					}
    				}
    				if (i == basisSize) {
    					basisIterator.setBasis(null);
    					super.reset();
    					return;
    				}
    			}
    			if (direction != null && !ignoreDirection &&
    					(agents[0].node.index() > agents[1].node.index()) == (direction == IteratorDirective.UP)) {
    				basisIterator.setBasis(null);
    				super.reset();
    				return;
    			}
    		}
    		basisIterator.setBasis(agents);
        	basisIterator.setTarget(targetAtoms);
    	}
    	super.reset();
    }
    
    public void setDirection(IteratorDirective.Direction direction) {
    	
    	if (!ignoreDirection) {//if this is true, basisIterator is instanceof ApiInterGroup
    		this.direction = direction;
    		needBasisUpdate = true;
    	}
    	else if (basisSize == 2) {
    		// if ignoreDirection is true and there are 2 species,
    		// basisIterator must be an intragroup iterator and
    		// the direction must be passed to it.  All other types
    		// of iterators do not need the direction.
    		((ApiIntragroup)basisIterator).setDirection(direction);
    	}
    }
    
    /**
     * Method used by constructor to determine the type of iterator to
     * wrap by this class.
     */
     private static AtomsetIterator makeIterator(Species[] species) {
    	if (species == null || species.length == 0 || species.length > 2
    			|| species[0] == null || species[species.length-1] == null) {
    		throw new IllegalArgumentException("null or invalid number of species.  Must specify either 1 or 2 species instances.");
    	}
    	if (species.length==1) {
    		return new AtomIteratorBasis();
    	}
    	else if (species[0] == species[1]) {
    		return new ApiIntragroup();
    	}
    	else return new ApiIntergroup();
    }
     
 	private final Species[] species;
	private final AtomsetIteratorBasisDependent basisIterator;
	private final SpeciesAgent[] agents;
	private final int basisSize;
	private boolean needBasisUpdate;
	private Atom[] targetAtoms;
	private IteratorDirective.Direction direction;
	private final boolean ignoreDirection;
	
}
