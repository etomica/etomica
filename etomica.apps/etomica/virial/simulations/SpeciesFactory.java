package etomica.virial.simulations;

import etomica.Simulation;
import etomica.Species;

public interface SpeciesFactory {

    public Species makeSpecies(Simulation sim);
}
