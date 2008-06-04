package etomica.virial;

import etomica.api.ISimulation;
import etomica.api.ISpecies;
import etomica.models.rowley.ConformationMethanol;
import etomica.models.rowley.SpeciesMethanol;
import etomica.space.ISpace;

/**
 * SpeciesFactory that makes methanol.
 */
public class SpeciesFactoryMethanol implements SpeciesFactory, java.io.Serializable {
	
	public SpeciesFactoryMethanol(boolean pointCharges) {
		this.pointCharges = pointCharges;
		
	}
    
    public ISpecies makeSpecies(ISimulation sim, ISpace space) {
    	
    	// The satellite site, X, is closer to the oxygen atom in the model with point charges.
    	SpeciesMethanol species = new SpeciesMethanol(space, pointCharges);
        
        return species;
    }
    
    private static final long serialVersionUID = 1L;
    protected boolean pointCharges;

}

