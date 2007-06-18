package etomica.virial;

import etomica.models.water.SpeciesWater3P;
import etomica.simulation.ISimulation;
import etomica.species.Species;


/**
 * SpeciesFactory that makes SpeciesWater
 */
public class SpeciesFactoryWater3P implements SpeciesFactory, java.io.Serializable {
    public Species makeSpecies(ISimulation sim) {
        return new SpeciesWater3P(sim);
    }
}
