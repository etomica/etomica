package etomica;

/**
 * Iterates over the neighbors of a particular atom, as specified by 
 * the atom's neighborManager, which is held as an allatomAgent.
 */
 
 /* History of changes
  * 09/16/02 (DAK) modified reset(Atom) method so that iteration is performed with
  *          current basis, beginning from given atom.  Previously this method would
  *          reset the basis to the given atom.
  */
public class AtomIteratorNeighbor implements AtomIterator {
    
    private NeighborManager neighborManager;
    private IteratorDirective.Direction direction = IteratorDirective.BOTH;
    private final AtomIteratorList iterator = new AtomIteratorList();
    public final int agentIndex;
    
    public AtomIteratorNeighbor() {
        this(Atom.requestAgentIndex(new Atom.AgentSource() {
            public Object makeAgent(Atom a) {
                return new NeighborManager(a);
            }
        }));
    }
    public AtomIteratorNeighbor(int agentIndex) {
        this.agentIndex = agentIndex;
 //       iterator.setSkipFirstAtom(true);//comment or not depending on whether Tab is used in NeighborManager
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
     * Sets iteration to begin with the given atom, using the current basis.
     * If given atom is not a neighbor of basis atom, hasNext is false.
     */
    public Atom reset(Atom atom) {
        /*deleted 09/16/02
        setBasis(atom);
        return iterator.reset(direction);
      added this:   */
        iterator.reset(neighborManager.tab, direction);
        //loop until given atom is found, or iterator expires
        while(iterator.hasNext() && iterator.peek() != atom) {iterator.next();}
        return iterator.peek();
    }

    public Atom reset(IteratorDirective.Direction direction) {
        return iterator.reset(neighborManager.tab, direction);
    }
    
    public void unset() {iterator.unset();}
    
    public Atom first() {
        throw new RuntimeException("method first() not implemented in AtomIteratorNeighbor");
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
        throw new RuntimeException("method AtomIteratorNeighbor.getBasis() not yet implemented");
    }
    public boolean contains(Atom a) {
        throw new RuntimeException("method AtomIteratorNeighbor.contains(Atom) not yet implemented");
    }
}//end of AtomIteratorNeighbor