package etomica.virial;

import etomica.api.ISpecies;
import etomica.space.ISpace;

public interface SpeciesFactory {

    public ISpecies makeSpecies(ISpace space);
}
