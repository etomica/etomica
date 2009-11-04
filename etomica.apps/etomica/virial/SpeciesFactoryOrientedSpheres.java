package etomica.virial;

import etomica.api.ISpecies;
import etomica.chem.elements.ElementSimple;
import etomica.space.ISpace;
import etomica.species.SpeciesSpheresRotating;


/**
 * SpeciesFactory that makes SpeciesSpheresMono 
 */
public class SpeciesFactoryOrientedSpheres implements SpeciesFactory, java.io.Serializable {
    public ISpecies makeSpecies(ISpace space) {
        return new SpeciesSpheresRotating(space, new ElementSimple("O"));
    }
}
