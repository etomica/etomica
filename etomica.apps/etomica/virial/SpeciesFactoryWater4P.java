package etomica.virial;

import etomica.config.Conformation;
import etomica.models.water.SpeciesWater4P;
import etomica.api.ISimulation;
import etomica.api.ISpecies;
import etomica.space.Space;


/**
 * SpeciesFactory that makes SpeciesWater
 */
public class SpeciesFactoryWater4P implements SpeciesFactory, java.io.Serializable {
    public SpeciesFactoryWater4P(Conformation conformation) {
        this.conformation = conformation;
    }
    public ISpecies makeSpecies(ISimulation sim, Space space) {
        SpeciesWater4P species = new SpeciesWater4P(space);
        species.getMoleculeType().setConformation(conformation);
        return species;
    }
    
    protected Conformation conformation;
}
