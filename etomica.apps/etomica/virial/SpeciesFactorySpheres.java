package etomica.virial;

import etomica.api.ISpecies;
import etomica.chem.elements.ElementSimple;
import etomica.space.ISpace;
import etomica.species.SpeciesSpheresMono;


/**
 * SpeciesFactory that makes SpeciesSpheresMono 
 */
public class SpeciesFactorySpheres implements SpeciesFactory, java.io.Serializable {
    public ISpecies makeSpecies(ISpace _space) {
        return new SpeciesSpheresMono(_space, new ElementSimple("A"));
    }
}
