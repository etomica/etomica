package etomica.potential;

import etomica.Atom;
import etomica.AtomsetIterator;
import etomica.IteratorDirective;
import etomica.Phase;
import etomica.Potential;
import etomica.Simulation;
import etomica.Space;
import etomica.IteratorDirective.Direction;
import etomica.atom.iterator.AtomsetIteratorBasisDependent;
import etomica.atom.iterator.AtomsetIteratorDirectable;
import etomica.atom.iterator.AtomsetIteratorTargetable;

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
     * Space is that for current value of Simulation.instance.
     */
    public PotentialGroup(int nBody) {
        this(nBody, Simulation.getDefault().space);
    }
    /**
     * Makes instance with null truncation, regardless of Default.TRUNCATE_POTENTIALS.
     */
    public PotentialGroup(int nBody, Space space) {
        this(nBody, space, null);
    }
    /**
     * Makes instance with given truncation scheme.
     */
    public PotentialGroup(int nBody, Space space, PotentialTruncation truncation) {
        super(nBody, space, truncation);
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
	protected synchronized void addPotential(Potential potential, AtomsetIterator iterator) {
		//the order of the given potential should be consistent with the order of the iterator
		if(potential.nBody() != iterator.nBody()) {
			throw new RuntimeException("Error: adding to PotentialGroup a potential and iterator that are incompatible");
		}
		//the given iterator should expect a basis of atoms equal in number to the order of this potential
		if(this.nBody() != ((AtomsetIteratorBasisDependent)iterator).basisSize()) {
			throw new RuntimeException("Error: adding an iterator that requires a basis size different from the nBody of the containing potential");
		}
		//put new potentials at beginning of list
		first = new PotentialLinker(potential, iterator, first);
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
	
	public double getRange() {
		return Double.MAX_VALUE;
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
				return;
			}//end if
			previous = link;
		}//end for
	}//end removePotential
    
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
    	iterator.reset();//loop over atom groups affected by this potential group
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
	
	protected PotentialLinker first;
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
    