package etomica;

/**
 * Iterates over the neighbors of a particular atom, as specified by 
 * the atom's neighborManager, which is held as an allatomAgent.
 */
 
 /* History of changes
  * 09/16/02 (DAK) modified reset(Atom) method so that iteration is performed with
  *          current basis, beginning from given atom.  Previously this method would
  *          reset the basis to the given atom.
  * 06/29/03 (DAK) added a flag (setBasisOnResetAtom) so that iterator can be
  * configured to interpret reset(Atom) as in modification above, or as the new
  * basis.
  */
public class AtomIteratorNeighbor implements AtomIterator {
    
    private NeighborManager neighborManager;
    private IteratorDirective.Direction direction = IteratorDirective.BOTH;
    private final AtomIteratorList iterator = new AtomIteratorList();
    public final int agentIndex;
    private final boolean setBasisOnResetAtom;
    
    /**
     * Constructs with new request for agent index from Atom, and
     * setBasisOnResetAtom = false.
     * @see java.lang.Object#Object()
     */
    public AtomIteratorNeighbor() {
        this(Atom.requestAgentIndex(new Atom.AgentSource() {
            public Object makeAgent(Atom a) {
                return new NeighborManager(a);
            }
        }),false);
    }
    /**
     * Constructs with given agent index, and setBasisOnResetAtom = false.
     * @param agentIndex see full-parameter constructor
     */
    public AtomIteratorNeighbor(int agentIndex) {
    	this(agentIndex, false);
    }
    
    /**
     * Constructs with new request for agent index from Atom, and given value of
     * setBasisOnResetAtom.
     * @param setBasisOnResetAtom see full-parameter constructor
     */
    public AtomIteratorNeighbor(boolean setBasisOnResetAtom) {
		this(Atom.requestAgentIndex(new Atom.AgentSource() {
			public Object makeAgent(Atom a) {
				return new NeighborManager(a);
			}
		}),setBasisOnResetAtom);
    }
    
    /**
     * Full-parameter constructor.
     * 
     * @param agentIndex indicates the index in the agents array of each atom
     * where neighbor manager for atom will be found
     * @param setBasisOnResetAtom flag affecting behavior of reset(Atom) and
     * reset(IteratorDirective) methods.  If true the given atom will be used to
     * set the basis, so that the iterator will return neighbors of the given
     * atom.  If false, uses reset method to indicate where looping is to begin,
     * taking the given atom to be a neighbor of the (already set) basis atom.
     * Default is false.  Set to true if using iterator as the inner-loop
     * iterator of a pair iterator.
     */
    public AtomIteratorNeighbor(int agentIndex, boolean setBasisOnResetAtom) {
        this.agentIndex = agentIndex;
        this.setBasisOnResetAtom = setBasisOnResetAtom;
 //       iterator.setSkipFirstAtom(true);//comment or not depending on whether Tab is used in NeighborManager
    }
    
	public void all(Atom basis, IteratorDirective id, final AtomActive action) {
		if(basis == null || basis.node.isLeaf() || action == null) return;
		iterator.all(basis, id, action);
//		throw new etomica.exception.MethodNotImplementedException();
	}

   public boolean hasNext() {return iterator.hasNext();}
    
    public Atom reset(IteratorDirective id) {
        direction = id.direction();
        switch(id.atomCount()) {
            case 0: return reset(direction);
            case 1: return reset(id.atom1());
            default: throw new IllegalArgumentException("AtomIteratorNeighbor.reset(IteratorDirective) unexpected atomCount");
        }
    }
    
    public Atom reset() {
        return reset(IteratorDirective.BOTH);
    }

    /**
     * If constructed with setBasisOnReset = false, then sets iteration to begin
     * with the given atom, using the current basis; otherwise, sets the given
     * atom as the basis and resets to iterate in current direction from it. If
     * setBasisOnReset = false and given atom is not a neighbor of basis atom,
     * hasNext is false.
     */
    public Atom reset(Atom atom) {
    	if(setBasisOnResetAtom) {
	        setBasis(atom);
	        return iterator.reset(((NeighborManager)atom.allatomAgents[agentIndex]).tab, direction);
    	} else {
	        iterator.reset(neighborManager.tab, direction);
			//loop until given atom is found, or iterator expires
		   while(iterator.hasNext() && iterator.peek() != atom) {iterator.next();}
		   return iterator.peek();
	   	}
    }

    public Atom reset(IteratorDirective.Direction direction) {
        return iterator.reset(neighborManager.tab, direction);
    }
    
    public void unset() {iterator.unset();}
    
    public Atom first() {
        throw new etomica.exception.MethodNotImplementedException();
    }
    public Atom next() {return iterator.next();}

    public void allAtoms(AtomAction act) {
        iterator.allAtoms(act);
    }
    
    public int size() {return neighborManager.neighborCount();}    
    
    public void setBasis(NeighborManager manager) {
        neighborManager = manager;
        iterator.setBasis(neighborManager.neighbors());
    }
    
    public void setBasis(Atom atom) {
        setBasis((NeighborManager)atom.allatomAgents[agentIndex]);
    }
    
    public void setupNeighbors(AtomList list, NeighborManager.Criterion criterion) {
     //   iterator.setBasis(list);
     //   iterator.reset();
        AtomIteratorListSimple iter = new AtomIteratorListSimple(list);
        while(iter.hasNext()) {
            Atom a = iter.next();
            NeighborManager manager = (NeighborManager)a.allatomAgents[agentIndex];
            manager.setupNeighbors(list, criterion);
        }
    }
    
    public Atom getBasis() {
        throw new etomica.exception.MethodNotImplementedException();
    }
    public boolean contains(Atom a) {
		throw new etomica.exception.MethodNotImplementedException();
    }
    
	/**
	 * Invokes all(Atom, IteratorDirective, AtomActive) method of this
	 * class, using given arguments if they are instances of the appropriate
	 * classes. Otherwise returns without throwing any exception.
	 * @see etomica.AtomSetIterator#all(AtomSet, IteratorDirective, AtomSetActive)
	 */
	public void all(AtomSet basis, IteratorDirective id, final AtomSetActive action) {
		 if(!(basis instanceof Atom && action instanceof AtomActive)) return;
		 all((Atom)basis, id, (AtomActive)action);
	}
    
    public static void main(String[] args) {
    	Simulation sim = Simulation.instance;
    	Species species = new SpeciesSpheresMono();
    	species.setNMolecules(5);
    	final Phase phase = new Phase();
		AtomIteratorNeighbor nbrIterator = new AtomIteratorNeighbor(true);
    	sim.elementCoordinator.go();
		NeighborManager.Criterion criterion = new NeighborManager.Criterion() {
			public boolean areNeighbors(Atom a1, Atom a2) {
				return Math.abs(a1.node.index()-a2.node.index()) == 1 ||
				   (a1==phase.firstAtom() && a2==phase.lastAtom()) ||
				   (a2==phase.firstAtom() && a1==phase.lastAtom());
			}};
		Atom first = phase.getAgent(species).node.firstLeafAtom();
		Atom middle = ((AtomTreeNodeGroup)phase.getAgent(species).node).childList.get(3);
		Atom last = phase.getAgent(species).node.lastLeafAtom();
		AtomList list = ((AtomTreeNodeGroup)phase.getAgent(species).node).childList;
		AtomIteratorList iterator = new AtomIteratorList(list);
		nbrIterator.setupNeighbors(list, criterion);
		nbrIterator.setBasis(first);
		IteratorDirective.testSuite(nbrIterator,first,middle,last);
		ApiGeneral api = new ApiGeneral(sim.space, iterator, nbrIterator);
		IteratorDirective.testSuitePair(api, first, middle, last);
    }
}//end of AtomIteratorNeighbor