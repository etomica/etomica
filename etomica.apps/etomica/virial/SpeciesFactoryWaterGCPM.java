package etomica.virial;

import etomica.models.water.ConformationWaterGCPM;
import etomica.models.water.SpeciesWater4P;
import etomica.api.ISimulation;
import etomica.api.ISpecies;
import etomica.space.Space;


/**
 * SpeciesFactory that makes SpeciesWater
 */
public class SpeciesFactoryWaterGCPM implements SpeciesFactory, java.io.Serializable {
    public ISpecies makeSpecies(ISimulation sim, Space _space) {
        SpeciesWater4P species = new SpeciesWater4P(_space);
        species.getMoleculeType().setConformation(new ConformationWaterGCPM(_space));
        return species;
    }
}
