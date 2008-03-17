package etomica.virial;

import etomica.api.ISimulation;
import etomica.api.ISpecies;
import etomica.space.Space;

public interface SpeciesFactory {

    public ISpecies makeSpecies(ISimulation sim, Space space);
}
