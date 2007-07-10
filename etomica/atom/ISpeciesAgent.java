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
     * Returns the Dimension (Quantity) corresponding to the NMolecules
     * property.
     */
    public Dimension getNMoleculesDimension();
}