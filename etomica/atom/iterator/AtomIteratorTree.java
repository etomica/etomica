package etomica.atom.iterator;

import etomica.action.AtomsetAction;
import etomica.action.AtomsetCount;
import etomica.atom.Atom;
import etomica.atom.AtomGroup;
import etomica.atom.AtomSet;

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
   
public class AtomIteratorTree implements AtomIterator, java.io.Serializable {
    
	/**
	 * Default gives a leaf-atom iterator.  Must set a root node and
	 * reset before using.
	 */
    public AtomIteratorTree() {
    	this(Integer.MAX_VALUE);
    }

    /**
     * Constructs iterator that will iterate over atoms at the given depth below
     * a (to-be-specified) root node.  Iterates atoms only at the given level, and
     * not those above it (doAllNodes = false by default).  Must set a root node
     * and reset before using.
     * @param d depth in tree for iteration.
     */
	public AtomIteratorTree(int d) {
	    this(null, d, false);
	}
    
    /**
     * Constructor permitting specification of all conditions.  Requires
     * reset before beginning iteration.
     * @param root          iteration root 
     * @param depth         nominal depth of iteration
     * @param doAllNodes    flag for iteration of all nodes between root and depth, inclusive
     */
    public AtomIteratorTree(Atom root, int depth, boolean doAllNodes) {
        setIterationDepth(depth);
        setRootAtom(root);
        setDoAllNodes(doAllNodes);
        counter = new AtomsetCount();
    }
    
    /**
     * Performs action on all iterates in current condition.  Unaffected
     * by reset status. Clobbers iteration state.
     */
    public void allAtoms(AtomsetAction act) {
        if(rootAtom == null) return;
		if(doAllNodes || iterationDepth == 0 || rootAtom.isLeaf()) {
            act.actionPerformed(rootAtom);
            if (iterationDepth == 0) return;
        }
		listIterator.reset();
		while(listIterator.hasNext()) {
            Atom atom = listIterator.nextAtom();
            if (atom.isLeaf() || iterationDepth == 1) {
                act.actionPerformed(atom);
            }
            else {
                if (treeIterator == null) {
                    treeIterator = new AtomIteratorTree(iterationDepth-1);
                    treeIterator.setDoAllNodes(doAllNodes);
                }
				treeIterator.setRootAtom(atom);
				treeIterator.allAtoms(act);
			}
		}
    	unset();
    }

    /**
     * Indicates if the iterator has another atom to give.
     */
    public boolean hasNext() {return next != null;}
    
    
    /**
     * Puts iterator in state in which hasNext is false.
     */
    public void unset() {next = null;}
    
    /**
     * Returns true if the given atom is among the iterates as
     * currently configured (independent of reset state).  Clobbers iteration state.
     */
    /**
     * Reinitializes the iterator according to the most recently specified basis,
     * iteration depth, and doAllNodes flag.
     */
    public void reset() {
        if(rootAtom == null) {
            unset();
            return;
        }
        listIterator.reset();
        if (treeIterator != null) treeIterator.unset();
        next = rootAtom;
        if(!doAllNodes && iterationDepth>0 && !rootAtom.isLeaf()) nextAtom();
    }

    /**
     * Returns the next atom in the iteration sequence.
     */
    public Atom nextAtom() {
        if(next == null) return null;
        Atom nextAtom = next;
        next = null;
        if (treeIterator != null && treeIterator.hasNext()) {
            next = treeIterator.nextAtom();
            return nextAtom;
        }
        while(listIterator.hasNext()) {
            Atom atom = listIterator.nextAtom();
            if (atom.isLeaf() || iterationDepth == 1) {
                next = atom;
                break;
            }
            if (treeIterator == null) {
                treeIterator = new AtomIteratorTree(iterationDepth-1);
                treeIterator.setDoAllNodes(doAllNodes);
            }
            treeIterator.setRootAtom(atom); 
            treeIterator.reset();
            if(treeIterator.hasNext()) {
                next = treeIterator.nextAtom();
                break;
            }
        }
        return nextAtom;
    }

    /**
     * Returns the next atom in the iteration sequence.  Same as nextAtom().
     */
    public AtomSet next() {
        return nextAtom();
    }

    /**
     * Returns the next atom without advancing the iterator.
     */
    public AtomSet peek() {
        return next;
    }
        
    /**
     * Defines the root of the tree under which the iteration is performed.
     * If atom is null, will return no atoms on reset, until a non-null root is specified.
     * User must perform a subsequent call to reset() before beginning iteration.
     */
    public void setRootAtom(Atom newRootAtom) {
        rootAtom = newRootAtom;
        if(newRootAtom != null) {
            if(iterationDepth == 0 || newRootAtom.isLeaf()) {//singlet iteration of basis atom
                if (!wealreadyknowyourstupid) {
                    System.err.println("don't use AtomIteratorTree as a singlet iterator.");
                    wealreadyknowyourstupid = true;
                }
                listIterator.setList(null);
            } else {
                listIterator.setList(((AtomGroup)rootAtom).getChildList());
            }
        }
        unset();
    }
        
    /**
     * Returns the number of iterates given by a full cycle of this iterator
     * in its current condition (independent of current iteration state).
     */
    public int size() {
        if(rootAtom == null) return 0;
    	unset();
        counter.reset();
        allAtoms(counter);
    	return counter.callCount();
    }
    
    /**
     * Returns 1, indicating that this is an Atom iterator.
     */
    public final int nBody() {return 1;}
    
    /**
     * Sets the depth below the current root for which the iteration will occur.
     * Any non-negative value is permitted.  A value of zero causes singlet iteration
     * returning just the root atom. A value of 1 returns all children of the root
     * atom, a value of 2 returns all children of all children of the root, etc.
     * If doAllNodes is false, iterator returns only atoms at the specified level, and not those above it.
     * Returns atoms at bottom of hierarchy (i.e., leafs) if specified depth
     * exceeds depth of hierarchy in a particular branch.  
     * Default is Integer.MAX_VALUE, which causes all leaf atoms to be iterated.
     */
    public void setIterationDepth(int depth) {
        if (iterationDepth == depth) return;
        if(depth < 0) throw new IllegalArgumentException("Error: attempt to set iteration depth to negative value in AtomIteratorTree");
        iterationDepth = depth;
        if(treeIterator != null && depth > 1) {
            treeIterator.setIterationDepth(depth - 1);
        }
        if(rootAtom != null) setRootAtom(rootAtom);
        unset();
    }
    /**
     * Returns the currently set value of iteration depth.
     */
    public int getIterationDepth() {return iterationDepth;}
    
    /**
     * Convenience method that sets iterator to iterate over all leaf atoms 
     * below the root atom.
     * Equivalent to setDoAllNodes(false) and setIterationDepth(Integer.MAX_VALUE)
     */
    public void setAsLeafIterator() {
        setDoAllNodes(false);
        setIterationDepth(Integer.MAX_VALUE);
    }
    
    /**
     * Accessor method for doAllNodes flag.
     */
	public boolean isDoAllNodes() {
		return doAllNodes;
	}
    
    /**
     * Flag indicating whether iterates are taken only at the iteration depth (false),
     * or whether all atoms between (and including) the root and iteration depth
     * are given (true).
     */
	public void setDoAllNodes(boolean doAllNodes) {
		this.doAllNodes = doAllNodes;
		if(treeIterator != null) treeIterator.setDoAllNodes(doAllNodes);
		unset();
	}

    private static final long serialVersionUID = 1L;
    private Atom rootAtom;
    private AtomIteratorTree treeIterator;//used for recursive iteration to lower levels in tree
    private final AtomIteratorArrayListSimple listIterator = new AtomIteratorArrayListSimple();
    private int iterationDepth = Integer.MAX_VALUE;
    private Atom next;
    private boolean doAllNodes = false;
    private boolean wealreadyknowyourstupid = false;
    private final AtomsetCount counter;
        
}//end of AtomIteratorTree
