package etomica.virial;

import etomica.simulation.ISimulation;
import etomica.species.ISpecies;
import etomica.species.SpeciesSpheresRotating;


/**
 * SpeciesFactory that makes SpeciesSpheresMono 
 */
public class SpeciesFactoryOrientedSpheres implements SpeciesFactory, java.io.Serializable {
    public ISpecies makeSpecies(ISimulation sim) {
        return new SpeciesSpheresRotating(sim);
    }
}
