package etomica.virial;

import etomica.models.water.SpeciesWater4P;
import etomica.api.IConformation;
import etomica.api.ISimulation;
import etomica.api.ISpecies;
import etomica.space.ISpace;


/**
 * SpeciesFactory that makes SpeciesWater
 */
public class SpeciesFactoryWater4P implements SpeciesFactory, java.io.Serializable {
    public SpeciesFactoryWater4P(IConformation conformation) {
        this.conformation = conformation;
    }
    public ISpecies makeSpecies(ISimulation sim, ISpace space) {
        SpeciesWater4P species = new SpeciesWater4P(space);
        species.setConformation(conformation);
        return species;
    }
    
    protected IConformation conformation;
}
