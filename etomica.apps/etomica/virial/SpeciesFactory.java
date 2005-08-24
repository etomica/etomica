package etomica.virial;

import etomica.Simulation;
import etomica.species.Species;

public interface SpeciesFactory {

    public Species makeSpecies(Simulation sim);
}
