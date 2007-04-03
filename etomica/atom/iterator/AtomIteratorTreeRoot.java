package etomica.atom.iterator;

import etomica.atom.Atom;

/**
 * Atom iterator that traverses all atoms at or to a specified depth below a
 * specified root atom in the atom tree hierarchy. The iterator is conditioned
 * using the following parameters:
 * <ul>
 * <li>root: atom that provides the point of reference for the iteration.  Atoms
 * below the iteration-root atom are subject to iteration.  The iteration-root
 * does not have to be the root of the species tree; any atom in the species tree
 * can serve as the root of iteration.
 * <li>depth: any non-negative integer, indicating how many levels below the root
 * can iterates be taken.  The value 0 indicates iteration of the root atom only, 
 * 1 is children of the root atom, etc. If the bottom of the hierarchy is reached 
 * before the specified depth, the leaf atoms encountered there are the iterates.
 * If the tree has branches that are deeper in some parts than in others, the deepest
 * iterates of each branch are taken up to the specified depth.  A very large value
 * of depth causes all leaf atoms below the root to be iterated.
 * <li>doAllNodes: flag indicating whether iteration is inclusive of all atoms
 * between the root and those at the iteration depth. If false, the atoms only 
 * at the specified depth (or leaf atoms if depth exceeds depth of tree) are iterated;
 * if true, all atoms from (and including) the root atom to those at the iteration
 * depth are iterated.  For example, if depth is 1, then doAllNodes == true indicates
 * iteration of both the root and its child atoms; doAllNodes == false indicates
 * only the child atoms are iterated. 
 * <ul>
 * 
 * @author David Kofke and Andrew Schultz
 */
   
public class AtomIteratorTreeRoot extends AtomIteratorTree {
    
	/**
	 * Default gives a leaf-atom iterator.  Must set a root node and
	 * reset before using.
	 */
    public AtomIteratorTreeRoot() {
    	this(Integer.MAX_VALUE);
    }

    /**
     * Constructs iterator that will iterate over atoms at the given depth below
     * a (to-be-specified) root node.  Iterates atoms only at the given level, and
     * not those above it (doAllNodes = false by default).  Must set a root node
     * and reset before using.
     * @param d depth in tree for iteration.
     */
	public AtomIteratorTreeRoot(int d) {
	    this(null, d, false);
	}
    
    /**
     * Constructor permitting specification of all conditions.  Requires
     * reset before beginning iteration.
     * @param root          iteration root 
     * @param depth         nominal depth of iteration
     * @param doAllNodes    flag for iteration of all nodes between root and depth, inclusive
     */
    public AtomIteratorTreeRoot(Atom root, int depth, boolean doAllNodes) {
        super(root, depth, doAllNodes);
    }
    
    public void setRootAtom(Atom newRootAtom) {
        // just make superclass method visible
        super.setRootAtom(newRootAtom);
    }

    private static final long serialVersionUID = 1L;
}
