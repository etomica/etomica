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
        if(nBody < 0 && !(this instanceof PotentialMaster)) throw new IllegalArgumentException("Cannot instantiate negative-body potential");
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
	 * Adds the given potential to this group, but should not be called directly.  Instead,
	 * this method is invoked by the setParentPotential method (or more likely, 
	 * in the constructor) of the given potential.  
	 */
    
    //this method needs given iterator to implement AtomsetDirectable, but PotentialMaster doesn't.
    //need to reconcile differences
	synchronized void addPotential(Potential potential, AtomsetIterator iterator) {
		if(potential == null || iterator == null) {
			throw new NullPointerException(); 
		}
		//the order of the given potential should be consistent with the order of the iterator
		if(potential.nBody() != iterator.nBody()) {
			throw new RuntimeException("Error: adding to PotentialGroup a potential and iterator that are incompatible");
		}
		//the given iterator should expect a basis of atoms equal in number to the order of this potential
/*		if(this.nBody() != ((AtomsetIteratorBasisDependent)iterator).basisSize()) {
			throw new RuntimeException("Error: adding an iterator that requires a basis size different from the nBody of the containing potential");
		}*/
		//Set up to evaluate zero-body potentials last, since they may need other potentials
		//to be configured for calculation (i.e., iterators set up) first
		if(((potential instanceof Potential0) || (potential instanceof PotentialGroupLrc)) && last != null) {//put zero-body potential at end of list
			last.next = makeLinker(potential, iterator, null);
			last = last.next;
		} else {//put other potentials at beginning of list
			first = makeLinker(potential, iterator, first);
			if(last == null) last = first;
		}
		addPotentialNotify(potential);
	}
	
	//TODO this needs some work
	public double energy(Atom[] basisAtoms) {
		if(basisAtoms.length != this.nBody()) {
			throw new IllegalArgumentException("Error: number of atoms for energy calculation inconsistent with order of potential");
		}
		double sum = 0.0;
		for (PotentialLinker link=first; link!= null; link=link.next) {		
			//if(firstIterate) ((AtomsetIteratorBasisDependent)link.iterator).setDirective(id);
			((AtomsetIteratorBasisDependent)link.iterator).setBasis(basisAtoms);
			link.iterator.reset();
			while(link.iterator.hasNext()) {
				sum += link.potential.energy(link.iterator.next());
			}
		}
		return sum;
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
		removePotentialNotify(potential);
	}//end removePotential
    
	public void addPotentialNotify(Potential potential) {
		if(parentElement instanceof PotentialGroup) {
			((PotentialGroup)parentElement).addPotentialNotify(potential);
		}
		else if(parentElement instanceof PotentialMaster) {//in case we make PotentialMaster not extend PotentialGroup
				((PotentialMaster)parentElement).addPotentialNotify(potential);
		}
		//else no notification is appropriate
	}
	
	
	public void removePotentialNotify(Potential potential) {
		if(parentElement instanceof PotentialGroup) {
			((PotentialGroup)parentElement).removePotentialNotify(potential);
		}
		else if(parentElement instanceof PotentialMaster) {//in case we make PotentialMaster not extend PotentialGroup
				((PotentialMaster)parentElement).removePotentialNotify(potential);
		}
		//else no notification is appropriate
	}
    /**
     * Performs the specified calculation over the iterates given by the iterator,
     * using the directive to set up the iterators for the sub-potentials of this group.
     */
	//TODO consider what to do with sub-potentials after target atoms are reached
    public void calculate(AtomsetIterator iterator, IteratorDirective id, PotentialCalculation pc) {
    	if(!enabled) return;
    	Atom[] targetAtoms = id.getTargetAtoms();
    	IteratorDirective.Direction direction = id.direction();
		//loop over sub-potentials
    	//TODO consider separate loops for targetable and directable
 		for (PotentialLinker link=first; link!= null; link=link.next) {
			((AtomsetIteratorTargetable)link.iterator).setTarget(targetAtoms);
			((AtomsetIteratorDirectable)link.iterator).setDirection(direction);
		}
    	iterator.reset();
		while (iterator.hasNext()) {
    		Atom[] basisAtoms = iterator.next();
    		for (PotentialLinker link=first; link!= null; link=link.next) {
    			((AtomsetIteratorBasisDependent)link.iterator).setBasis(basisAtoms);   			
    			pc.doCalculation(link.iterator, id, link.potential);
    		}
     	}
    }//end calculate
    
    public void setPhase(Phase phase) {
    	this.phase = phase;
  		for (PotentialLinker link=first; link!= null; link=link.next) {
			link.potential.setPhase(phase);
		}
    }
    
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
	protected boolean enabled = true;
	protected Phase phase;

	protected static class PotentialLinker implements java.io.Serializable {
	    protected final Potential potential;
	    protected final AtomsetIterator iterator;
	    protected PotentialLinker next;
	    //Constructors
	    public PotentialLinker(Potential a, AtomsetIterator i, PotentialLinker l) {
	    	potential = a;
	    	iterator = i;
	    	next = l;
	    }
	}

}//end PotentialGroup
    