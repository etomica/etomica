package etomica.atom;

import etomica.api.IAtom;
import etomica.api.IAtomSet;
import etomica.api.IAtomType;
import etomica.api.IMolecule;

public class Molecule extends Atom implements IMolecule {

    public Molecule(IAtomType type) {
        super(type);
        childList = new AtomArrayList();
        assignChildOrdinals();
    }
    
    /**
     * Returns a string of digits that uniquely identifies this atom.  String is
     * formed by concatenating the ordinal of this atom to the signature
     * given by the parent of this atom.  If atom has no parent, forms a string
     * from only the ordinal.
     */
    public String signature() {
        return Integer.toString(index);
    }

    /**
     * Returns a string formed by concatenating the signature of this atom
     * to a string that identifies it as a molecule.
     */
    public final String toString() {
//        return Integer.toBinaryString(node.index());
        return "Molecule(" + index + ")";
    }

    /**
     * Assigns ordinals to all child atoms, numbering them sequentially
     * according to their position in the childList.
     */
    private void assignChildOrdinals() {
        for (int i = 0; i < childList.getAtomCount(); i++) {
            childList.getAtom(i).setIndex(i);
        }
    }

    /**
     * Adds the given Atom as a child of this Atom.  The given child Atom
     * should be parentless when this method is called.
     * @throws IllegalArgumentException if the given atom already has a parent.
     */
    public void addChildAtom(IAtomLeaf newChildAtom) {
        if(newChildAtom.getParentGroup() != null) {//new parent is null
            throw new IllegalArgumentException(newChildAtom+" is already the child of "+newChildAtom.getParentGroup());
        }

        newChildAtom.setParent(this);

        newChildAtom.setIndex(childList.getAtomCount());
        childList.add(newChildAtom);
    }
    
    /**
     * Removes the given child Atom from this AtomGroup.
     * @throws IllegalArgumentException if the given atom is not a child.
     */
    public void removeChildAtom(IAtomLeaf oldChildAtom) {
        for (int i=0; i<childList.getAtomCount(); i++) {
            if (childList.getAtom(i) == oldChildAtom) {
                oldChildAtom.setParent(null);
                childList.removeAndReplace(i);
                childList.maybeTrimToSize();
                if (childList.getAtomCount() > i) {
                    // reassign the old last Atom (which is now in the removed
                    // Atom's place) to have the old Atom's index.
                    childList.getAtom(i).setIndex(i);
                }
                return;
            }
        }
        throw new IllegalArgumentException(oldChildAtom+" is not a child");
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
    public IAtom getDescendant(int[] path) {
        return getDescendant(0, path);
    }
    
    private IAtom getDescendant(int n, int[] path) {
        IAtom child = childList.getAtom(path[n]);
        if(path.length - 1 > n) {//go further down hierarchy
            if(!(child instanceof IMolecule)) {//no more there
                throw new IllegalArgumentException("Depth of requested descendant exceeds depth of atom hierarchy");
            }//get indicated descendant recursively
            child = ((Molecule)child).getDescendant(n+1, path);
        }
        return child;
    }
    
    /**
     * @return the childList
     */
    public final IAtomSet getChildList() {
        return childList;
    }

    private static final long serialVersionUID = 1L;
    
    //nobody should add/remove atoms except via AtomGroup's methods.
    //consider a mechanism to ensure this; a inner mutator class made available only
    //to list's creator, for example (still wouldn't prevent modification via direct
    //access of entry classes).
    protected final AtomArrayList childList;
}
