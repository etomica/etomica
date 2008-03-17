package etomica.virial;

import etomica.api.ISimulation;
import etomica.api.ISpecies;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;


/**
 * SpeciesFactory that makes SpeciesSpheresMono 
 */
public class SpeciesFactorySpheres implements SpeciesFactory, java.io.Serializable {
    public ISpecies makeSpecies(ISimulation sim, Space _space) {
        return new SpeciesSpheresMono(sim, _space);
    }
}
