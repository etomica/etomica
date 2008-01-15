package etomica.virial;

import etomica.models.water.ConformationWaterGCPM;
import etomica.models.water.SpeciesWater4P;
import etomica.simulation.ISimulation;
import etomica.species.Species;


/**
 * SpeciesFactory that makes SpeciesWater
 */
public class SpeciesFactoryWaterGCPM implements SpeciesFactory, java.io.Serializable {
    public Species makeSpecies(ISimulation sim) {
        SpeciesWater4P species = new SpeciesWater4P(sim);
        species.getMoleculeType().setConformation(new ConformationWaterGCPM(sim.getSpace()));
        return species;
    }
}
