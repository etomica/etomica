package etomica.virial;

import etomica.models.water.SpeciesWater3P;
import etomica.simulation.Simulation;
import etomica.species.Species;


/**
 * SpeciesFactory that makes SpeciesWater
 */
public class SpeciesFactoryWater3P implements SpeciesFactory, java.io.Serializable {
    public Species makeSpecies(Simulation sim) {
        return new SpeciesWater3P(sim);
    }
}
