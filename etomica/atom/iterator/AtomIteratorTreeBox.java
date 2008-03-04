package etomica.atom.iterator;

import etomica.action.AtomAction;
import etomica.action.AtomsetAction;
import etomica.api.IAtom;
import etomica.api.IAtomSet;
import etomica.api.IBox;
import etomica.api.IMolecule;

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
   
public class AtomIteratorTreeBox extends AtomIteratorTree implements AtomIteratorBoxDependent {
    
	/**
	 * Default gives a leaf-atom iterator.
	 */
    public AtomIteratorTreeBox() {
    	this(Integer.MAX_VALUE);
    }

    /**
     * Constructs iterator that will iterate over atoms at the given depth below
     * a (to-be-specified) root node.  Iterates atoms only at the given level, and
     * not those above it (doAllNodes = false by default).
     * @param d depth in tree for iteration.
     */
	public AtomIteratorTreeBox(int d) {
	    this(null, d, false);
	}
    
    /**
     * Constructor permitting specification of all conditions.  Requires
     * reset before beginning iteration.
     * @param box         box containing atoms to iterator over 
     * @param depth         nominal depth of iteration
     * @param doAllNodes    flag for iteration of all nodes between root and depth, inclusive
     */
    public AtomIteratorTreeBox(IBox box, int depth, boolean doAllNodes) {
        super(null, depth, doAllNodes);
        if (box != null) {
            setBox(box);
        }
    }
    
    public void setIterationDepth(int newDepth) {
        if (newDepth < 1) {
            throw new RuntimeException("depth must be at least 1");
        }
        super.setIterationDepth(newDepth);
    }
    
    public void setBox(IBox newBox) {
        box = newBox;
        listIterator.setList(newBox.getMoleculeList());
        unset();
    }
    
    /**
     * Performs action on all iterates in current condition.  Unaffected
     * by reset status. Clobbers iteration state.
     */
    public void allAtoms(AtomsetAction act) {
        listIterator.setList(box.getMoleculeList());
        listIterator.reset();
        for (IAtom atom = listIterator.nextAtom(); atom != null;
             atom = listIterator.nextAtom()) {
            if (!(atom instanceof IMolecule) || iterationDepth == 1) {
                atomSetSinglet.atom = atom;
                act.actionPerformed(atomSetSinglet);
                continue;
            }
            
            if (doAllNodes) {
                atomSetSinglet.atom = atom;
                act.actionPerformed(atomSetSinglet);
            }

            IAtomSet childList = ((IMolecule)atom).getChildList();
            for (int iChild=0; iChild<childList.getAtomCount(); iChild++) {
                atomSetSinglet.atom = childList.getAtom(iChild);
                act.actionPerformed(atomSetSinglet);
            }
        }
        unset();
    }

    /**
     * Performs action on all iterates in current condition.  Unaffected
     * by reset status. Clobbers iteration state.
     */
    public void allAtoms(AtomAction act) {
        listIterator.setList(box.getMoleculeList());
        listIterator.reset();
        for (IAtom atom = listIterator.nextAtom(); atom != null;
             atom = listIterator.nextAtom()) {
            if (!(atom instanceof IMolecule) || iterationDepth == 1) {
                act.actionPerformed(atom);
                continue;
            }
            
            if (doAllNodes) {
                act.actionPerformed(atom);
            }

            IAtomSet childList = ((IMolecule)atom).getChildList();
            for (int iChild=0; iChild<childList.getAtomCount(); iChild++) {
                act.actionPerformed(childList.getAtom(iChild));
            }
        }
        unset();
    }

    public void reset() {
        listIterator.setList(box.getMoleculeList());
        listIterator.reset();

        if (childListIterator != null) childListIterator.unset();
    }

    /**
     * Returns the number of iterates given by a full cycle of this iterator
     * in its current condition (independent of current iteration state).
     */
    public int size() {
        if(box == null) return 0;
        unset();
        counter.reset();
        allAtoms(counter);
        return counter.callCount();
    }

    private static final long serialVersionUID = 1L;
    protected IBox box;
}
