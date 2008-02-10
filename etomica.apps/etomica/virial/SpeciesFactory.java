package etomica.virial;

import etomica.simulation.ISimulation;
import etomica.species.ISpecies;

public interface SpeciesFactory {

    public ISpecies makeSpecies(ISimulation sim);
}
