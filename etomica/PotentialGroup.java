package etomica;

/**
 * Collection of potentials that act between the atoms contained in
 * two groups of atoms.  This group iterates over all such atom-group
 * pairs assigned to it.  For each pair it iterates over the potentials it
 * contains, instructing these sub-potentials to perform their calculations
 * over the atoms relevant to them in the two groups.
 *
 * @author David Kofke
 */

 /* History of changes
  * 8/14/02 (DAK) introduced removePotential method and modified addPotential.
  */

public class PotentialGroup extends Potential {
    
    /**
     * Constructor for use only by PotentialMaster subclass.
	 * @param sim Simulation instance in which potential is used.
     */
    PotentialGroup(Simulation sim) {
    	super(sim);
    }
    /**
     * Makes instance with null truncation, regardless of Default.TRUNCATE_POTENTIALS.
     * Parent potential is potential master for current value of Simulation.instance.
     */
    public PotentialGroup() {
        this(Simulation.instance.hamiltonian.potential);
    }
    /**
     * Makes instance with null truncation, regardless of Default.TRUNCATE_POTENTIALS.
     */
    public PotentialGroup(PotentialGroup parent) {
        this(parent, null);
    }
    /**
     * Makes instance with given truncation scheme.
     */
    public PotentialGroup(PotentialGroup parent, PotentialTruncation truncation) {
        super(parent, truncation);
    }

    public boolean contains(Potential potential) {
        for(PotentialLinker link=first; link!=null; link=link.next) {
            if(link.potential.equals(potential)) return true;
        }//end for
        return false;
    }
 
	/**
	 * Adds the given potential to this group, but should not be called directly.  Instead,
	 * this method is invoked by the setParentPotential method (or more likely, 
	 * in the constructor) of the given potential.  
	 */
	synchronized void addPotential(Potential potential) {
		if(potential == null || potential.parentPotential() != this) 
			throw new IllegalArgumentException("Improper call to addPotential; should use setParentPotential method in child instead of addPotential method in parent group"); 
		//Set up to evaluate zero-body potentials last, since they may need other potentials
		//to be configured for calculation (i.e., iterators set up) first
		if((potential instanceof Potential0) || (potential instanceof PotentialGroupLrc) && last != null) {//put zero-body potential at end of list
			last.next = makeLinker(potential, null);
			last = last.next;
		} else {//put other potentials at beginning of list
			first = makeLinker(potential, first);
			if(last == null) last = first;
		}
	}
	
	/**
	 * Returns a linker that is used to form linked list of potentials.  May be
	 * overridden to permit use of specialized linkers that hold additional
	 * information about the potential (this is done by SpeciesMaster).
	 * @param p the potential
	 * @param next the linker that is to follow the new one
	 * @return PotentialLinker the new potential linker
	 */
	protected PotentialLinker makeLinker(Potential p, PotentialLinker next) {
		return new PotentialLinker(p, next);
	}

	
	/**
	 * Removes given potential from the group.  No error is generated if
	 * potential is not in group.
	 */
	public synchronized void removePotential(Potential potential) {
		PotentialLinker previous = null;
		for(PotentialLinker link=first; link!=null; link=link.next) {
			if(link.potential == potential) {//found it
				if(previous == null) first = link.next;  //it's the first one
				else previous.next = link.next;          //it's not the first one
				if(link == last) last = previous; //removing last; this works also if last was also first (then removing only, and set last to null)
				return;
			}//end if
			previous = link;
		}//end for
	}//end removePotential
	   
	PotentialLinker first, last;
    
    /**
     * Performs the specified calculation over the iterates of this potential
     * that comply with the iterator directive.
     */
    public void calculate(AtomSet basis, IteratorDirective id, PotentialCalculation pc) {
    	if(!enabled) return;
		PotentialCalculation.PotentialGroupWrapper wrapper = pc.wrapper().set(this);
		wrapper.localDirective.copy(id);
		iterator.all(basis, wrapper.localDirective, wrapper);
		wrapper.release();
    }//end calculate
    
    AtomSetIterator iterator;
    
}//end PotentialGroup
    