package etomica.virial;

import etomica.models.water.SpeciesWater4P;
import etomica.simulation.ISimulation;
import etomica.species.Species;


/**
 * SpeciesFactory that makes SpeciesWater
 */
public class SpeciesFactoryWater4P implements SpeciesFactory, java.io.Serializable {
    public Species makeSpecies(ISimulation sim) {
        return new SpeciesWater4P(sim);
    }
}
