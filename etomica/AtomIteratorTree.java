package etomica;

/**
 * Atom iterator that traverses the elements of a tree of atoms.
 */
 
public class AtomIteratorTree implements AtomIterator {
    
    public boolean hasNext() {return next != null;}
    
    public boolean contains(Atom atom) {
        throw new RuntimeException("contains method not implemented in AtomIteratorTree");
    }
    
    public Atom reset(IteratorDirective id) {
        throw new RuntimeException("allAtoms method not implemented in AtomIteratorTree");
    }
    
    public Atom reset() {
        listIterator.reset();
        if(doTreeIteration) 
            do {
                treeIterator.reset(listIterator.next());
            } while(!treeIterator.hasNext() && listIterator.hasNext());
        return next = iterator.next();
    }
    
    public Atom reset(Atom atom) {
        if(atom == null) return next = null;
        listIterator.setBasis(((AtomTreeNodeGroup)atom.node).childList);
        return reset();
    }
    
    /**
     * Returns the next atom in the iteration sequence.
     */
    public Atom next() {
        Atom nextAtom = next;
        if(iterator.hasNext()) next = iterator.next();
 //       else if(listIterator.hasNext()) {//if this is evaluated and is true, then iterator == treeIterator
        else {
            do {
                if(!listIterator.hasNext()) {
                    next = null;
                    return nextAtom;
                } //else iterator == treeIterator
                if(basisIsMaster) {
                    treeIterator.setBasis(listIterator.next());
                    treeIterator.reset();
                } else {//not iterating over Agents -- assume structure of all molecules is the same
                    treeIterator.reset(listIterator.next());
                }
            } while(!treeIterator.hasNext());
            next = treeIterator.next();
        }
        return nextAtom;
    }
    
    public void allAtoms(AtomAction act) {
        throw new RuntimeException("allAtoms method not implemented in AtomIteratorTree");
    }
    
    /**
     * Defines generally the atoms subject to iteration.
     */
    public void setBasis(Atom atom) {
        basis = (AtomTreeNodeGroup)atom.node;
        AtomList list = basis.childList;
        basisIsMaster = (atom instanceof SpeciesMaster);
        listIterator.setBasis(list);
        doTreeIteration = (iterationDepth > 1 && list.size() > 0 && basis.childrenAreGroups());
        if(doTreeIteration) {
            //if doing tree iteration
            if(treeIterator == null) treeIterator = new AtomIteratorTree();
            treeIterator.setIterationDepth(iterationDepth-1);
            treeIterator.setBasis(list.getFirst());
            iterator = treeIterator;
        } else {
            iterator = listIterator;
        }
    }
    
    public Atom getBasis() {return basis.atom();}
        
    public int size() {
        throw new RuntimeException("size method not implemented in AtomIteratorTree");
    }
    
    public void setIterationDepth(int depth) {iterationDepth = depth;}
    public int getIterationDepth() {return iterationDepth;}
    
    public void setAsNeighbor(boolean b) {
        throw new RuntimeException("setAsNeighbor not implemented in AtomIteratorTree");
    }
    
    private AtomTreeNodeGroup basis;
    private int iteratorDepth;
    private boolean doTreeIteration, basisIsMaster;
    private AtomIteratorTree treeIterator;
    private final AtomIteratorList listIterator = new AtomIteratorList();
    private AtomIterator iterator;
    private boolean hasNext;
    private int iterationDepth = Integer.MAX_VALUE;
    private Atom next;
    
    public static void main(String args[]) {
        
        Simulation sim = new Simulation();
        Simulation.instance = sim;
        Species species2 = new SpeciesSpheresMono();
        Species species1 = new SpeciesSpheres(3,3);
        Species species0 = new SpeciesSpheres(3,2);
        species0.setNMolecules(0);
        species1.setNMolecules(2);
        species2.setNMolecules(3);
        Phase phase = new Phase();
        sim.elementCoordinator.go();
        
        AtomIteratorTree treeIterator = new AtomIteratorTree();
 //       treeIterator.setIterationDepth(2);
        treeIterator.setBasis(phase.speciesMaster);
        treeIterator.reset();
        while(treeIterator.hasNext()) {
            System.out.println(treeIterator.next().toString());
        }
    }
    
}