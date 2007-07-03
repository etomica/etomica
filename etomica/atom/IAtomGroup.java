package etomica.atom;

/**
 * Interface for a group of IAtoms, typically a molecule or a SpeciesAgent.
 */
public interface IAtomGroup extends IAtom {

    /**
     * Adds the given Atom as a child of this Atom.  The given child Atom
     * should be parentless when this method is called.
     * @throws IllegalArgumentException if the given atom already has a parent.
     */
    public void addChildAtom(IAtom newChildAtom);

    /**
     * Removes the given child Atom from this AtomGroup.
     * @throws IllegalArgumentException if the given atom is not a child.
     */
    public void removeChildAtom(IAtom oldChildAtom);

    /**
     * Notifies this atom group that the given child IAtom has been added to it
     * or one of its descendants.  This method should be called by addChildAtom
     * or by this IAtomGroup's child if the child is an IAtomGroup and has a
     * child added to it.
     */
    public void addAtomNotify(IAtom childAtom);

    /**
     * Notifies this atom group that the given child IAtom has been removed
     * from it or one of its descendants.  This method should be called by
     * removeChildAtom or by this IAtomGroup's child if the child is an
     * IAtomGroup and has a child removed from it.
     */
    public void removeAtomNotify(IAtom childAtom);

    /**
     * @return the children as an AtomArrayList
     */
    public AtomSet getChildList();
}