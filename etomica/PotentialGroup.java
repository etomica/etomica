package etomica;

/**
 * Collection of potentials that act between the atoms contained in
 * one or more groups of atoms.  This group iterates over all such atom-groups
 * assigned to it.  For each group it iterates over the potentials it contains,
 * instructing these sub-potentials to perform their calculations over the atoms
 * relevant to them in the groups.
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
    public PotentialGroup(int nBody) {
        this(nBody, Simulation.instance.hamiltonian.potential);
    }
    /**
     * Makes instance with null truncation, regardless of Default.TRUNCATE_POTENTIALS.
     */
    public PotentialGroup(int nBody, SimulationElement parent) {
        this(nBody, parent, null);
    }
    /**
     * Makes instance with given truncation scheme.
     */
    public PotentialGroup(int nBody, SimulationElement parent, PotentialTruncation truncation) {
        super(nBody, parent, truncation);
        switch(nBody) {
        	case 1: iterator = simulation().iteratorFactory.makeGroupIteratorSequential();
        			break;
        	case 2: iterator = new Api1A(new ApiIntergroup1A(simulation()),
											new ApiIntergroupAA(simulation()));
					break;
			default: throw new etomica.exception.MethodNotImplementedException("Potential group not developed to handle more than one- or two-body interactions");
        }//end switch
    }

	/**
	 * Indicates if a given potential is a sub-potential of this group.
	 * @param potential the potential in question
	 * @return boolean true if potential has been added to this group
	 */
    public boolean contains(Potential potential) {
        for(PotentialLinker link=first; link!=null; link=link.next) {
            if(link.potential.equals(potential)) return true;
        }//end for
        return false;
    }

	/**
	  * Extends superclass setSpecies method to instantiate iterators based on
	  * number of species in array.  If this is a one-body potential, iterator
	  * is instantiated as in Potential1, if a two-body, iterator is as in
	  * Potential2.
	  * @see etomica.Potential#setSpecies(Species[])
	  * @see etomica.Potential1#setSpecies(Species[])
	  * @see etomica.Potential2#setSpecies(Species[])
	  */
	 public void setSpecies(Species[] species) {
		 super.setSpecies(species);
		 switch(nBody) {
		 	case 1:
		 			break;
		 	case 2:
				 if(species.length == 1 || species[0] == species[1]) {
					 iterator = new Api1A(new ApiIntragroup1A(simulation()),
											 new ApiIntragroupAA(simulation()));
				 } else {
					 iterator = new Api1A(new ApiIntergroup1A(simulation()),
											 new ApiIntergroupAA(simulation()));
				 }
				 break;
		 }//end switch
	 }
 
	/**
	 * Adds the given potential to this group, but should not be called directly.  Instead,
	 * this method is invoked by the setParentPotential method (or more likely, 
	 * in the constructor) of the given potential.  
	 */
	synchronized void addPotential(Potential potential, AtomsetIterator iterator) {
		if(potential == null || potential.parentPotential() != this) 
			throw new IllegalArgumentException("Improper call to addPotential; should use setParentPotential method in child instead of addPotential method in parent group"); 
		//Set up to evaluate zero-body potentials last, since they may need other potentials
		//to be configured for calculation (i.e., iterators set up) first
		if(((potential instanceof Potential0) || (potential instanceof PotentialGroupLrc)) && last != null) {//put zero-body potential at end of list
			if (potential instanceof PotentialGroup) {
				lastGroup.next = makeLinker(potential, iterator, null);
				lastGroup = lastGroup.next;
			}
			else {
				last.next = makeLinker(potential, null);
				last = last.next;
			}
		} else {//put other potentials at beginning of list
			if (potential instanceof PotentialGroup) {
				firstGroup = makeLinker(potential, iterator, firstGroup);
				if (lastGroup == null) lastGroup = firstGroup;
			}
			else {
				first = makeLinker(potential, iterator, first);
				if(last == null) last = first;
			}
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
	protected PotentialLinker makeLinker(Potential p, AtomsetIterator i, PotentialLinker next) {
		return new PotentialLinker(p, i, next);
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
    
    /**
     * Performs the specified calculation over the iterates of this potential
     * that comply with the iterator directive.
     */
    //wrapper is used to avoid problems that might arise if more than one thread
    //is accessing the potential.  Wrapper is obtained from the potential calculation,
    //which assumes that it would be exclusive to a given thread.  Wrapper is
    //a PotentialCalculation subclass with an actionPerformed method that loops
    //over the sub-potentials of this PotentialGroup.
    public void calculate(AtomsetIterator iterator, IteratorDirective id, PotentialCalculation pc) {
    	if(!enabled) return;
		//loop over sub-potentials
    	iterator.reset();
    	while (iterator.hasNext()) {
    		Atom[] basisAtoms = iterator.next();
    		for (PotentialLinker link=first; link!= null; link=link.next) {
    			link.iterator.setBasis(basisAtoms);
    			link.iterator.setDirective(id);
    			pc.doCalculation(link.iterator,link.potential);
    		}
    		for (PotentialLinker link=firstGroup; link!=null; link=link.next) {
    			link.iterator.setBasis(basisAtoms);
    			link.iterator.setDirective(id);
    			link.potential.calculate(basisAtoms, id, pc);
    		}
    	}
    }//end calculate
    
    
	/**
	 * Returns the enabled flag, which if false will cause potential to not
	 * contribute to any potential calculations.
	 * @return boolean The current value of the enabled flag
	 */
	public boolean isEnabled() {
		return enabled;
	}

	/**
	 * Sets the enabled flag, which if false will cause potential to not
	 * contribute to any potential calculations.
	 * @param enabled The enabled value to set
	 */
	public void setEnabled(boolean enabled) {
		this.enabled = enabled;
	}
	
	protected PotentialLinker first, last;
	protected PotentialLinker firstGroup, lastGroup;
	protected boolean enabled = true;

	protected static class PotentialLinker implements java.io.Serializable {
	    protected final Potential potential;
	    protected final AtomsetIterator iterator;
	    protected PotentialLinker next;
	    //Constructors
//	    public PotentialLinker(Potential a) {potential = a;}
	    public PotentialLinker(Potential a, AtomsetIterator i, PotentialLinker l) {
	    	potential = a;
	    	iterator = i;
	    	next = l;
	    }
	}

}//end PotentialGroup
    