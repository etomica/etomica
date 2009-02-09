package etomica.api;

import etomica.chem.elements.Element;

public interface IAtomTypeLeaf extends IAtomType {

    public void setSpecies(ISpecies newParent);

    /**
     * Informs the ISpecies what its index should be.  This should only be
     * called by the SpeciesManager.
     */
    public void setIndex(int newIndex);

    /**
     * Returns the index for this IAtomTypeLeaf, within the context of an
     * ISimulation.  The index is the IAtomTypeLeaf' position in the list of
     * species and atom types in the simulation.
     */
    public int getIndex();
    
    public void setChildIndex(int newChildIndex);

    public int getChildIndex();

    public ISpecies getSpecies();

    /**
     * Returns the value of the mass.
     */
    public double getMass();

    /**
     * Returns the reciprocal of the mass, 1.0/mass
     */
    public double rm();

    public Element getElement();

}