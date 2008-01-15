package etomica.virial;

import etomica.config.Conformation;
import etomica.models.water.SpeciesWater4P;
import etomica.simulation.ISimulation;
import etomica.species.Species;


/**
 * SpeciesFactory that makes SpeciesWater
 */
public class SpeciesFactoryWater4P implements SpeciesFactory, java.io.Serializable {
    public SpeciesFactoryWater4P(Conformation conformation) {
        this.conformation = conformation;
    }
    public Species makeSpecies(ISimulation sim) {
        SpeciesWater4P species = new SpeciesWater4P(sim);
        species.getMoleculeType().setConformation(conformation);
        return species;
    }
    
    protected Conformation conformation;
}
