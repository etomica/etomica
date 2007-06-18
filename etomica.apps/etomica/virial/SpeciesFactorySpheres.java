package etomica.virial;

import etomica.simulation.ISimulation;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;


/**
 * SpeciesFactory that makes SpeciesSpheresMono 
 */
public class SpeciesFactorySpheres implements SpeciesFactory, java.io.Serializable {
    public Species makeSpecies(ISimulation sim) {
        return new SpeciesSpheresMono(sim);
    }
}
