package etomica.virial;

import etomica.simulation.ISimulation;
import etomica.species.Species;

public interface SpeciesFactory {

    public Species makeSpecies(ISimulation sim);
}
