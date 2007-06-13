package etomica.atom;

import etomica.units.Dimension;

public interface ISpeciesAgent extends IAtomGroup {

    /**
     * Returns the AtomManager responsible for managing this SpeciesAgent
     */
    public AtomManager getAtomManager();

    /**
     * Returns the number of molecules associated with this SpeciesAgent (this
     * will be equal to the number of IAtoms in the childList).
     */
    public int getNMolecules();

    /**
     * Sets the number of molecules for this species.  Molecules are either
     * added or removed until the given number is obtained.  Takes no action
     * at all if the new number of molecules equals the existing number.
     *
     * @param n  the new number of molecules for this species
     */
    public void setNMolecules(int n);

    /**
     * Returns the Dimension (Quantity) corresponding to the NMolecules
     * property.
     */
    public Dimension getNMoleculesDimension();

    /**
     * Adds a new molecule to this SpeciesAgent (it will be added to the
     * childList).
     */
    public IAtom addNewAtom();
}