package etomica.virial;

import etomica.Simulation;
import etomica.Species;
import etomica.SpeciesSpheresMono;


/**
 * SpeciesFactory that makes SpeciesSpheresMono 
 */
public class SpeciesFactorySpheres implements SpeciesFactory, java.io.Serializable {
    public Species makeSpecies(Simulation sim) {
        return new SpeciesSpheresMono(sim);
    }
}
