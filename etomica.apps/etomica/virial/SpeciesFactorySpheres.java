package etomica.virial;

import etomica.api.ISimulation;
import etomica.api.ISpecies;
import etomica.species.SpeciesSpheresMono;


/**
 * SpeciesFactory that makes SpeciesSpheresMono 
 */
public class SpeciesFactorySpheres implements SpeciesFactory, java.io.Serializable {
    public ISpecies makeSpecies(ISimulation sim) {
        return new SpeciesSpheresMono(sim);
    }
}
