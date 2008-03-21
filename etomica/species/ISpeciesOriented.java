package etomica.species;

import etomica.api.ISpecies;
import etomica.api.IVector;

public interface ISpeciesOriented extends ISpecies {

    /**
     * Returns the principle components of the moment of inertia of the
     * molecule within the body-fixed frame.  Do NOT modify the returned moment
     * of inertia returned.
     */
    public IVector getMomentOfInertia();

    public double getMass();

}