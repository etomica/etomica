package etomica.virial;

import etomica.models.water.SpeciesWater3P;
import etomica.api.ISimulation;
import etomica.api.ISpecies;


/**
 * SpeciesFactory that makes SpeciesWater
 */
public class SpeciesFactoryWater3P implements SpeciesFactory, java.io.Serializable {
    public ISpecies makeSpecies(ISimulation sim) {
        return new SpeciesWater3P(sim.getSpace());
    }
}
