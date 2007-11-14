package etomica.virial;

import etomica.simulation.ISimulation;
import etomica.species.Species;
import etomica.species.SpeciesSpheresRotating;


/**
 * SpeciesFactory that makes SpeciesSpheresMono 
 */
public class SpeciesFactoryOrientedSpheres implements SpeciesFactory, java.io.Serializable {
    public Species makeSpecies(ISimulation sim) {
        return new SpeciesSpheresRotating(sim);
    }
}
