package etomica.atom;


/**
 * Implementation of AtomTreeNode for non-leaf node.
 */
  
public class AtomTreeNodeGroup extends AtomTreeNode {
    
    public AtomTreeNodeGroup(Atom atom) {
        super(atom);
        childList = new AtomArrayList();
        assignChildOrdinals();
    }
    
    /**
     * Set this atom's index and update indexes of its descendants.
     */
    protected void setIndex(int parentIndex, int ordinal) {
        super.setIndex(parentIndex, ordinal);
        if (childList != null) {
            assignChildOrdinals();
        }
    }

    /**
     * Assigns ordinals to all child atoms, numbering them sequentially
     * according to their position in the childList.
     */
    private void assignChildOrdinals() {
        for (int i = 0; i < childList.size(); i++) {
            childList.get(i).getNode().setIndex(i);
        }
    }

    /**
     * Returns a specified atom descended from this one in the atom tree.  
     * Each index of the given array specifies the i-th child at the
     * depth of the array index.  So if path is {2, 0, 3},
     * returns the 3rd child of the 0th child of the 2nd child of
     * this node.  That is: (this node) -> (2nd child) -> (0th child) -> (3rd child)
     * The path indexes do not correspond to the ordinals assigned to the
     * children (ordinals are numbered from 1; specifications in path are
     * numbered from 0).
     */
    public Atom getDescendant(int[] path) {
        return getDescendant(0, path);
    }
    
    private Atom getDescendant(int n, int[] path) {
        Atom child = childList.get(path[n]);
        if(path.length - 1 > n) {//go further down hierarchy
            if(child.getNode().isLeaf()) {//no more there
                throw new IllegalArgumentException("Depth of requested descendant exceeds depth of atom hierarchy");
            }//get indicated descendant recursively
            child = ((AtomTreeNodeGroup)child.getNode()).getDescendant(n+1, path);
        }
        return child;
    }
	
    public boolean isLeaf() {return false;}
    
    public AtomLeaf firstLeafAtom() {
        for (int i = 0; i < childList.size(); i++) {
            AtomLeaf a1 = childList.get(i).getNode().firstLeafAtom();
            if(a1 != null) return a1;
        }
        return null;
    }
    
    /**
     * Returns the last leaf atom descended from this group.
     */
    public AtomLeaf lastLeafAtom() {
        for (int i = childList.size()-1; i > -1; i--) {
            AtomLeaf a1 = childList.get(i).getNode().lastLeafAtom();
            if(a1 != null) return a1;
        }
        return null;
    }
    
    public int leafAtomCount() {return leafAtomCount;}
    public int childAtomCount() {return childList.size();}

    /**
     * Returns the children of this group in an array of atoms.
     * Array is constructed on-the-fly, and is not updated with any
     * subsequent atom addition/removals.  Since array construction is
     * involved, this method can be expensive in computationally intensive
     * situations involving repeated calls (this should be avoided).
     */
    public Atom[] childAtomArray() {
        return childList.toArray();
    }
    
    public void removeAllChildren() {
        Atom[] array = childAtomArray();
        for(int i=0; i<array.length; i++) {
            array[i].getNode().dispose();
        }
    }
    
    /**
     * Notifies this atom group that an atom has been added to it 
     * or one of its descendants.
     */
    public void addAtomNotify(Atom childAtom) {
        leafAtomCount += childAtom.getNode().leafAtomCount();
        if (parentNode() != null) {
            parentNode().addAtomNotify(childAtom);
        }
    }
    
    /**
     * Notifies this atom group that an atom has been removed from it or 
     * one of its descendants.
     */
    public void removeAtomNotify(Atom childAtom) {
        leafAtomCount -= childAtom.getNode().leafAtomCount();
        if(parentNode() != null) {
            parentNode().removeAtomNotify(childAtom);
        }
    }

    /**
     * @return the childList
     */
    public final AtomArrayList getChildList() {
        return childList;
    }

    private static final long serialVersionUID = 1L;
    protected int leafAtomCount;
    
    //childlist is public, but should not add/remove atoms except via node's methods.
    //consider a mechanism to ensure this; a inner mutator class made available only
    //to list's creator, for example (still wouldn't prevent modification via direct
    //access of entry classes).
    private final AtomArrayList childList;
        
    public static final AtomTreeNodeFactory FACTORY = new AtomTreeNodeGroup.Factory();
    
    protected static final class Factory implements AtomTreeNodeFactory, java.io.Serializable {
        public AtomTreeNode makeNode(Atom atom) {
            return new AtomTreeNodeGroup(atom);
        }
        private static final long serialVersionUID = 1L;
    }
}