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
	 * Part of the PotentialGroup interface.
	 */
	public synchronized void addPotential(Potential potential) {
		if(potential == null || potential.parentPotential() != this) 
			throw new IllegalArgumentException("Improper call to addPotential; should use setParentPotential method in child instead of addPotential method in parent group"); 
		first = new PotentialLinker(potential, first);
		agentIterator.reset();
		while(agentIterator.hasNext()) {
			((Agent)agentIterator.next()).addPotential(potential.requestAgent());
		}
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
			}
			previous = link;
		}
		agentIterator.reset();
		while(agentIterator.hasNext()) {
			((Agent)agentIterator.next()).removePotential(potential);
		}

	}
	
	/**
	 * Makes a new agent group and fills its list with new agents from all
	 * potentials presently in this potential group.
	 * @see etomica.Potential#makeAgent()
	 */
	public synchronized PotentialAgent makeAgent() {
		Agent agent = new Agent();
		for(PotentialLinker link=first; link!=null; link=link.next) {
			agent.addPotential(link.potential.makeAgent());
		}//end for
		return agent;		
	}
   
	private PotentialLinker first;

	private class Agent extends PotentialAgent implements AtomSetAction {
    	
		private final IteratorDirective localDirective = new IteratorDirective();
		private final PotentialAgent.List groupList = new PotentialAgent.List();
		private final PotentialAgent.Iterator groupIterator = groupList.iterator();
		private PotentialCalculation potentialCalculation;
		
		private Agent() {
			super(PotentialGroup.this);
		}
    
	    /**
	     * Performs the specified calculation over the iterates of this potential
	     * that comply with the iterator directive.
	     */
	    public void calculate(IteratorDirective id, PotentialCalculation pc) {
	        potentialCalculation = pc;
			localDirective.copy(id);//copy the iteratordirective to define the directive sent to the subpotentials
			iterator.all(basis, localDirective, this);
	    }//end calculate
	    
	    public void action(AtomSet atomSet) {
			if(potentialTruncation != null && potentialTruncation.isZero(atomSet)) return;                
	            
			//if the atom of the pair is the one specified for calculation, then
			//it becomes the basis for the sub-potential iterations, and is no longer
			//specified to them via the iterator directive
			if(atomSet.contains(id.atom1())) localDirective.set();
	            
			//loop over sub-potentials
			groupIterator.reset();
			while(groupIterator.hasNext()) {
				PotentialAgent potential = groupIterator.next();
				if(localDirective.excludes(potential)) continue; //see if potential is ok with iterator directive
				potential.set(atomSet).calculate(localDirective, potentialCalculation);
			}//end for	    	
	    }
	
	    /**
	     * Convenient reformulation of the calculate method, applicable if the potential calculation
	     * performs a sum.  The method returns the potential calculation object, so that the sum
	     * can be accessed in-line with the method call.
	     */
	    public final PotentialCalculation.Sum calculate(IteratorDirective id, PotentialCalculation.Sum pa) {
	        this.calculate(id, (PotentialCalculation)pa);
	        return pa;
	    }
	    
	    /**
	     * Adds the given potential to this group, but should not be called directly.  Instead,
	     * this method is invoked by the setParentPotential method (or more likely, 
	     * in the constructor) of the given potential.  
	     * Part of the PotentialGroup interface.
	     */
	    private void addPotential(PotentialAgent agent) {
	    	groupList.add(agent);	    	
	    }
	
	    /**
	     * Removes given potential from the group.  No error is generated if
	     * potential is not in group.
	     */
	    private void removePotential(Potential potential) {
	    	groupIterator.reset();
	    	while(groupIterator.hasNext()) {
	    		PotentialAgent agent = groupIterator.next();
	    		if(agent.potential == potential) groupList.remove(agent);
	    	}
	    }
    }//end of Agent
    
}//end PotentialGroup
    