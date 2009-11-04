package etomica.models.traPPE;

import etomica.api.ISpecies;
import etomica.space.ISpace;
import etomica.virial.SpeciesFactory;

/**
 * SpeciesFactory that makes methanol.
 */
public class SpeciesFactoryMethanol implements SpeciesFactory, java.io.Serializable {
	
	public SpeciesFactoryMethanol() {
		
	}
    
    public ISpecies makeSpecies(ISpace space) {
    	
    	// The satellite site, X, is closer to the oxygen atom in the model with point charges.
    	SpeciesMethanol species = new SpeciesMethanol(space);
        
        return species;
    }
    
    private static final long serialVersionUID = 1L;

}


