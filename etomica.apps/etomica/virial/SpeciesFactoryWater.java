package etomica.virial;

import etomica.Simulation;
import etomica.Species;
import etomica.models.water.SpeciesWater;


/**
 * SpeciesFactory that makes SpeciesWater
 */
public class SpeciesFactoryWater implements SpeciesFactory, java.io.Serializable {
    public Species makeSpecies(Simulation sim) {
        return new SpeciesWater(sim);
    }
}
