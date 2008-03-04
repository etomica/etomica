package etomica.virial;

import etomica.api.ISimulation;
import etomica.api.ISpecies;

public interface SpeciesFactory {

    public ISpecies makeSpecies(ISimulation sim);
}
