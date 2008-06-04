package etomica.virial;

import etomica.api.ISimulation;
import etomica.api.ISpecies;
import etomica.models.rowley.SpeciesMethanol;
import etomica.space.ISpace;

/**
 * SpeciesFactory that makes methanol.
 */
public class SpeciesFactoryMethanol implements SpeciesFactory, java.io.Serializable {
    
    public ISpecies makeSpecies(ISimulation sim, ISpace space) {
        SpeciesMethanol species = new SpeciesMethanol(space);
        return species;
    }
    
    private static final long serialVersionUID = 1L;

}

