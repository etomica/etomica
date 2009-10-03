package etomica.api;


/**
 * Interface for an IMolecule
 */
public interface IMolecule {

    /**
     * Returns this IMolecule's index, which is its place in the list of
     * molecules of its ISpecies in its IBox.
     */
    public int getIndex();

    /**
     * Informs the IMolecule of its index.
     */
    public void setIndex(int index);

    /**
     * Returns the atoms in the molecule as an IAtomList.
     */
    public IAtomList getChildList();

    /**
     * Returns the ISpecies of this IMolecule.
     */
    public ISpecies getType();
}