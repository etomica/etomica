package etomica.atom;

public interface IAtom extends AtomSet {

    /**
     * Returns true if this atom is in the same molecule as the given atom.
     * 
     * @throws NullPointerException
     *             if the argument is null
     */
    public boolean inSameMolecule(IAtom atom);

    public void setGlobalIndex(SpeciesMaster speciesMaster);

    public int getGlobalIndex();

    /**
     * @return the Atom type, holding properties held in common with other 
     * atoms made by this atom's factory.
     */
    public AtomType getType();

    public boolean isLeaf();

    /**
     * Informs the Atom that the given AtomGroup is its parent.
     * This method should only be called by the parent.
     */
    public void setParent(AtomGroup newParent);

    public AtomGroup getParentGroup();

    /**
     * Integer assigned to this atom by its parent molecule.
     */
    public int getAddress();

    public void setIndex(int index);

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