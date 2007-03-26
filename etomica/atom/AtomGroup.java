package etomica.atom;

import java.util.ArrayList;

public class AtomGroup extends Atom {

    public AtomGroup(AtomType type) {
        super(type);
        childList = new AtomArrayList();
        assignChildOrdinals();
    }
    
    /**
     * Set this atom's index and update indexes of its descendants.
     */
    public void setIndex(int parentIndex) {
        super.setIndex(parentIndex);
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
            childList.get(i).setIndex(i);
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
            if(child.isLeaf()) {//no more there
                throw new IllegalArgumentException("Depth of requested descendant exceeds depth of atom hierarchy");
            }//get indicated descendant recursively
            child = ((AtomGroup)child).getDescendant(n+1, path);
        }
        return child;
    }
    
    public boolean isLeaf() {return false;}
    
    public int childAtomCount() {return childList.size();}

    /**
     * Returns the children of this group in an array of atoms.
     * Array is constructed on-the-fly, and is not updated with any
     * subsequent atom addition/removals.  Since array construction is
     * involved, this method can be expensive in computationally intensive
     * situations involving repeated calls (this should be avoided).
     */
    public Atom[] childAtomArray() {
        ArrayList foo = new ArrayList();
        foo.toArray();
        return childList.toArray();
    }
    
    /**
     * Notifies this atom group that an atom has been added to it 
     * or one of its descendants.
     */
    public void addAtomNotify(Atom childAtom) {
        if (parent != null) {
            parent.addAtomNotify(childAtom);
        }
    }
    
    /**
     * Notifies this atom group that an atom has been removed from it or 
     * one of its descendants.
     */
    public void removeAtomNotify(Atom childAtom) {
        if(parent != null) {
            parent.removeAtomNotify(childAtom);
        }
    }

    /**
     * @return the childList
     */
    public final AtomArrayList getChildList() {
        return childList;
    }

    private static final long serialVersionUID = 1L;
    
    //nobody should not add/remove atoms except via AtomGroup's methods.
    //consider a mechanism to ensure this; a inner mutator class made available only
    //to list's creator, for example (still wouldn't prevent modification via direct
    //access of entry classes).
    protected final AtomArrayList childList;
}
