package etomica.virial;

import etomica.api.ISimulation;
import etomica.api.ISpecies;
import etomica.space.ISpace;

public interface SpeciesFactory {

    public ISpecies makeSpecies(ISimulation sim, ISpace space);
}
