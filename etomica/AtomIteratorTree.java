package etomica;

/**
 * Atom iterator that traverses all atoms at a specified depth below a
 * specified basis atom in the atom tree hierarchy.  The depth may be
 * specified as any non-negative integer.  If the bottom of the hierarchy
 * is reached before the specified depth, the leaf atoms encountered
 * there are the iterates.  Assumes that all molecules (i.e., atoms
 * below the SpeciesAgent level) of a given species have the same structure.
 * For example, there is not a species in which one molecule has two levels of
 * atoms below it, while another molecule of the same species has one level
 * of atoms below it.  However, does not assume that molecules in different species 
 * have the same structure.
 *
 * @author David Kofke
 * 02.02.16
 */
 
public class AtomIteratorTree implements AtomIterator {
    
    public AtomIteratorTree() {}
    
    /**
     * Returns new iterator ready to loop over all the leaf atoms below
     * the given atom.
     */
    public AtomIteratorTree(Atom atom) {
        this();
        setBasis(atom);
        reset();
    }
    
    /**
     * Returns a new iterator ready to loop over all the atoms at the
     * given depth below the given atom.
     */
    public AtomIteratorTree(Atom atom, int d) {
        this();
        setIterationDepth(d);
        setBasis(atom);
        reset();
    }
    
    /**
     * Indicates if the iterate has another atom to give.
     */
    public boolean hasNext() {return next != null;}
    
    /**
     * Not yet implemented.
     */
    public boolean contains(Atom atom) {
        throw new RuntimeException("contains method not implemented in AtomIteratorTree");
    }
    
    /**
     * Not implemented.
     */
    public Atom reset(IteratorDirective id) {
        throw new RuntimeException("allAtoms method not implemented in AtomIteratorTree");
    }
    
    /**
     * Reinitializes the iterator according to the most recently specified basis
     * and iteration depth.
     */
    public Atom reset() {
        listIterator.reset();
        if(doTreeIteration) 
            do {
                if(basisIsMaster) {
                    treeIterator.setBasis(listIterator.next()); 
                    treeIterator.reset();
                }
                else treeIterator.reset(listIterator.next());
            } while(!treeIterator.hasNext() && listIterator.hasNext());
        return next = iterator.next();
    }
    
    /**
     * Resets the iterator using the given atom as the basis.  Method is
     * designed for use only by another AtomIteratorTree class.  As a
     * tree iterator loops over children of the basis atom, it calls
     * this reset method of another tree iterator instance to perform
     * iterations over the tree atoms below the child.  This assumes that
     * each child atom in this sequence has the same tree structure below it.
     */
    Atom reset(Atom atom) {
        if(atom == null || basis == null) return next = null;
        listIterator.setBasis(((AtomTreeNodeGroup)atom.node).childList);
        return reset();
    }
    
    /**
     * Returns the next atom in the iteration sequence.
     */
    public Atom next() {
        Atom nextAtom = next;
        if(iterator.hasNext()) next = iterator.next();
        else {
            do {//advance list iterator and reset treeIterator until list expires or tree hasNext
                if(!listIterator.hasNext()) {
                    next = null;
                    return nextAtom;
                } //else iterator == treeIterator, because !iterator.hasNext() so iterator != listIterator
                if(basisIsMaster) {//looping over species agents, so use setBasis to prepare tree for different structure of molecules
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
    
    /**
     * Not yet implemented.
     */
    public void allAtoms(AtomAction act) {
        throw new RuntimeException("allAtoms method not implemented in AtomIteratorTree");
    }
    
    /**
     * Defines the root of the tree under which the iteration is performed.
     * If atom is null, hasNext will report false.
     * User must perform a subsequent call to reset() before beginning iteration.
     */
    public void setBasis(Atom atom) {
        if(atom == null) {
            basis = null; 
            listIterator.setBasis(AtomList.NULL); 
            doTreeIteration = false;
            return;
        }
        basis = (AtomTreeNodeGroup)atom.node;
        basisIsMaster = (atom instanceof SpeciesMaster);
        if(iterationDepth == 0 || atom.node.isLeaf()) {//singlet iteration of basis atom
            if(singletList == null) singletList = new AtomList();
            singletList.clear();
            singletList.add(atom);
            listIterator.setBasis(singletList);
            doTreeIteration = false;
            return;
        }
        AtomList list = basis.childList;
        listIterator.setBasis(list);
        doTreeIteration = (iterationDepth > 1 &&
                            (basisIsMaster || (list.size() > 0 && basis.childrenAreGroups())));
        if(doTreeIteration) {
            if(treeIterator == null) treeIterator = new AtomIteratorTree();
            treeIterator.setIterationDepth(iterationDepth-1);
            treeIterator.setBasis(list.getFirst());
            iterator = treeIterator;
        } else {
            iterator = listIterator;
        }
    }
    
    /**
     * Returns the current root of the iteration tree.
     */
    public Atom getBasis() {return basis.atom();}
    
    /**
     * Not yet implemented.
     */
    public int size() {
        throw new RuntimeException("size method not implemented in AtomIteratorTree");
    }
    
    /**
     * Sets the depth below the current basis for which the iteration will occur.
     * Any non-negative value is permitted.  A value of zero causes singlet iteration
     * returning just the basis atom. A value of 1 returns all children of the basis
     * atom, a value of 2 returns all children of all children of the basis, etc.
     * Iterator returns only atoms at the specified level, and not those above it.
     * Returns atoms at bottom of hierarchy (i.e., leaves) if specified depth exceeds 
     * depth of hierarchy.  Defaults is Integer.MAX_VALUE, which causes all leaf atoms 
     * to be iterated.
     */
    public void setIterationDepth(int depth) {
        if (iterationDepth == depth) return;
        if(depth < 0) throw new IllegalArgumentException("Error: attempt to set iteration depth to negative value in AtomIteratorTree");
        iterationDepth = depth;
        if(basis != null) setBasis(getBasis());
    }
    /**
     * Returns the currently set value of iteration depth.
     */
    public int getIterationDepth() {return iterationDepth;}
    
    /**
     * Sets iterator to iterate over all leaf atoms below the basis atom.
     * Equivalent to setIterationDepth(Integer.MAX_VALUE).
     */
    public void setLeafIterator() {
        setIterationDepth(Integer.MAX_VALUE);
    }
    
    /**
     * Not implemented.
     */
    public void setAsNeighbor(boolean b) {
        throw new RuntimeException("setAsNeighbor not implemented in AtomIteratorTree");
    }
    
    private AtomTreeNodeGroup basis;
    private int iteratorDepth;
    private boolean doTreeIteration, basisIsMaster;
    private AtomIteratorTree treeIterator;
    private final AtomIteratorList listIterator = new AtomIteratorList();
    private AtomIterator iterator;
    private AtomList singletList;
    private boolean hasNext;
    private int iterationDepth = Integer.MAX_VALUE;
    private Atom next;
    
    /**
     * main method to test and demonstrate use of this class.
     */
    public static void main(String args[]) {
        
        Simulation sim = new Simulation();
        Simulation.instance = sim;
        Species species2 = new SpeciesSpheresMono();
        Species species1 = new SpeciesSpheres(3,3);
        Species species0 = new SpeciesSpheres(3,2);
        species0.setNMolecules(0);
        species1.setNMolecules(2);
        species2.setNMolecules(2);
        Phase phase = new Phase();
        sim.elementCoordinator.go();
        
        AtomIteratorTree treeIterator = new AtomIteratorTree();
 //       treeIterator.setIterationDepth(2);
        treeIterator.setBasis(phase.speciesMaster);
        treeIterator.reset();
        while(treeIterator.hasNext()) System.out.println(treeIterator.next().toString());
        System.out.println();
        
        treeIterator.setIterationDepth(2);
        treeIterator.reset();
        while(treeIterator.hasNext()) System.out.println(treeIterator.next().toString());
        System.out.println();
        
        treeIterator.setIterationDepth(1);
        treeIterator.reset();
        while(treeIterator.hasNext()) System.out.println(treeIterator.next().toString());
        System.out.println();
        
        treeIterator.setIterationDepth(0);
        treeIterator.reset();
        while(treeIterator.hasNext()) System.out.println(treeIterator.next().toString());
        System.out.println();
        
        treeIterator.setBasis(phase.speciesMaster.node.firstChildAtom());
        treeIterator.setLeafIterator();
        treeIterator.setIterationDepth(1);
        treeIterator.reset();
        while(treeIterator.hasNext()) System.out.println(treeIterator.next().toString());
        System.out.println();
        
        System.out.print("null-basis test ");
        treeIterator.setBasis(null);
        treeIterator.setLeafIterator();
        treeIterator.reset();
        while(treeIterator.hasNext()) System.out.println(treeIterator.next().toString());
        System.out.println("ok");
    }//end of main
    
}//end of AtomIteratorTree