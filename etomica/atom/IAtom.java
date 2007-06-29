package etomica.atom;

public interface IAtom extends AtomSet {

    /**
     * Returns true if this atom is in the same molecule as the given atom.
     * 
     * @throws NullPointerException
     *             if the argument is null
     */
    public boolean inSameMolecule(IAtom atom);

    /**
     * Informs the IAtom that it should ask the SpeciesMaster for a global
     * index (with the Box).  This method should only be called by the
     * SpeciesMaster.
     */
    public void setGlobalIndex(AtomManager atomManager);

    /**
     * Returns the global index (within the Box) of this Atom.  The global
     * index is unique to the IAtom in the Box.  The IAtom's global may
     * change over the course of a simulation due to addition or removal of
     * other IAtoms in the Box.  An BoxGlobalAtomIndexEvent is fired by
     * the Box's event manager when an Atom's global index changes. 
     */
    public int getGlobalIndex();

    /**
     * @return the Atom type, holding properties held in common with other 
     * atoms made by this atom's factory.
     */
    public AtomType getType();

    /**
     * Informs the Atom that the given AtomGroup is its parent.
     * This method should only be called by the parent.
     */
    public void setParent(IAtomGroup newParent);

    /**
     * Returns the parent AtomGruop of this IAtom.
     */
    public IAtomGroup getParentGroup();

    /**
     * Returns the address indicating the position of this IAtom within the
     * Box's IAtom hierarchy.  The address is a bitmasked 32-byte integer
     * with a set of bits appropriated to each level of the IAtom hierarchy.
     * The address is typically interpreted by an AtomAddressManager.
     */
    public int getAddress();

    /**
     * Informs the IAtom of its index, which is used to construct the address.
     */
    public void setIndex(int index);

    /**
     * Returns this IAtom's index, which is its place in the parent AtomGroup's
     * list of child IAtoms.
     */
    public int getIndex();

    /**
     * Returns true if the given atom is in the hierarchy of parents of this atom,
     * or if the given atom is this atom.  Returns true, for example, if the given
     * atom is this atom's parent, or its parent's parent, etc.
     */
    public boolean isDescendedFrom(IAtom group);

    /**
     * Returns the child of the given node through which this node is derived.
     * If given node is parent node of this, returns this.
     * If this node is not descended from the given node, returns null.
     */
    public IAtom getChildWhereDescendedFrom(IAtom atom);

}