package etomica.virial;

import etomica.api.ISimulation;
import etomica.api.ISpecies;
import etomica.space.Space;
import etomica.species.SpeciesSpheresRotating;


/**
 * SpeciesFactory that makes SpeciesSpheresMono 
 */
public class SpeciesFactoryOrientedSpheres implements SpeciesFactory, java.io.Serializable {
    public ISpecies makeSpecies(ISimulation sim, Space space) {
        return new SpeciesSpheresRotating(sim, space);
    }
}
