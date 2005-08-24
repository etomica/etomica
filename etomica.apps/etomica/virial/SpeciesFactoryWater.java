package etomica.virial;

import etomica.models.water.SpeciesWater;
import etomica.simulation.Simulation;
import etomica.species.Species;


/**
 * SpeciesFactory that makes SpeciesWater
 */
public class SpeciesFactoryWater implements SpeciesFactory, java.io.Serializable {
    public Species makeSpecies(Simulation sim) {
        return new SpeciesWater(sim);
    }
}
