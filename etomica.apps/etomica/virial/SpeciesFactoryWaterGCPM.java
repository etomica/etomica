package etomica.virial;

import etomica.models.water.ConformationWaterGCPM;
import etomica.models.water.SpeciesWater4P;
import etomica.api.ISimulation;
import etomica.api.ISpecies;


/**
 * SpeciesFactory that makes SpeciesWater
 */
public class SpeciesFactoryWaterGCPM implements SpeciesFactory, java.io.Serializable {
    public ISpecies makeSpecies(ISimulation sim) {
        SpeciesWater4P species = new SpeciesWater4P(sim.getSpace());
        species.getMoleculeType().setConformation(new ConformationWaterGCPM(sim.getSpace()));
        return species;
    }
}
