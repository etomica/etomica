/*
 * History
 * Created on Aug 31, 2004 by kofke and schultz
 */
package etomica;


//this class is slated for deletion


/**
 * Iterates molecules (children of SpeciesAgent atoms), as specified by phase
 * and species.  Species must be specified at construction, and when a 
 * subsequent call to setPhase is made, the basis for iteration will be 
 * set to the species agents in that phase.  
 */
//PotentialMaster uses instances of this class for all iterators it assigns
//to the potentials added to it via setSpecies.
public final class AtomsetIteratorMoleculeOriginal extends AtomsetIteratorAdapter
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
	public AtomsetIteratorMoleculeOriginal(Species[] species) {
        this(species,IteratorFactorySimple.INSTANCE);
    }
    
    public AtomsetIteratorMoleculeOriginal(Species[] species, IteratorFactory iteratorFactory) {
		super(iteratorFactory.makeMoleculeIterator(species));
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
    	//TODO this is probably inefficient.  swatting a fly with a sledgehammer
    	System.arraycopy(targetAtoms,0,this.targetAtoms,0,Math.min(this.targetAtoms.length,targetAtoms.length));
    	for(int i=targetAtoms.length; i<this.targetAtoms.length; i++) {
    		this.targetAtoms[i] = null;
    	}
    	needBasisUpdate = true;
    	unset();
    }
    
    public void reset() {
        //target and/or phase changed and phase is not null
    	if (needBasisUpdate && agents[0] != null) {
        	needBasisUpdate = false;
            //target has been set
    		if (targetAtoms != null && targetAtoms[0] != null) {
    			SpeciesAgent targetAgent = targetAtoms[0].node.parentSpeciesAgent();
    			if (targetAgent != agents[0]) {
    				int i=basisSize;
    				if (!ignoreDirection) {//inter-species iteration
                        //look for target in other species
    					for (i=1; i<basisSize; i++) {
    						if (targetAgent == agents[i]) {//found it; make target agent 0
    							agents[i] = agents[0];
    							agents[0] = targetAgent;
    							break;
    						}
    					}
    				}
                    //ignoring direction, or target not related to this iterator's species
    				if (i == basisSize) {
    					basisIterator.setBasis(null);
    					super.reset();
    					return;
    				}
    			}
                //targetAgent is now agent[0]
                //set up iteration in relation to direction specification
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
    	else if (basisSize == 2 && basisIterator instanceof AtomsetIteratorDirectable) {
    		// if ignoreDirection is true and there are 2 species,
    		// basisIterator must be an intragroup iterator and
    		// the direction must be passed to it.  All other types
    		// of iterators do not need the direction.
    		((AtomsetIteratorDirectable)basisIterator).setDirection(direction);
    	}
    }
     
 	private final Species[] species;
	private final AtomsetIteratorBasisDependent basisIterator;
	private final SpeciesAgent[] agents;
	private final int basisSize;
	private boolean needBasisUpdate;
	private final Atom[] targetAtoms;
	private IteratorDirective.Direction direction;
	private final boolean ignoreDirection;
	
}
