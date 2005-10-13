package etomica.virial;

import etomica.models.water.SpeciesWater4P;
import etomica.simulation.Simulation;
import etomica.species.Species;


/**
 * SpeciesFactory that makes SpeciesWater
 */
public class SpeciesFactoryWater4P implements SpeciesFactory, java.io.Serializable {
    public Species makeSpecies(Simulation sim) {
        return new SpeciesWater4P(sim);
    }
}
