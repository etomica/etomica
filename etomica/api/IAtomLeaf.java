package etomica.api;


public interface IAtomLeaf extends IAtom {

    /**
     * Sets the atom's global index to the give value.  This method should only
     * be called by the IBox.
     */
    public void setLeafIndex(int newGlobalIndex);

    /**
     * Returns the global index (within the Box) of this Atom.  The global
     * index is unique to the IAtom in the Box.  The IAtom's global may
     * change over the course of a simulation due to addition or removal of
     * other IAtoms in the Box.  An BoxGlobalAtomIndexEvent is fired by
     * the Box's event manager when an Atom's global index changes. 
     */
    public int getLeafIndex();

    /**
     * Informs the Atom that the given AtomGroup is its parent.
     * This method should only be called by the parent.
     */
    public void setParent(IMolecule newParent);

    /**
     * Returns the parent AtomGruop of this IAtom.
     */
    public IMolecule getParentGroup();

}
