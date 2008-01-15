package etomica.atom;

public interface IAtomLeaf extends IAtom {

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
