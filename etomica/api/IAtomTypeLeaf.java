package etomica.api;

import etomica.chem.elements.Element;
import etomica.units.Dimension;

public interface IAtomTypeLeaf extends IAtomType {

    public void setSpecies(ISpecies newParent);

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